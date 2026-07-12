import boto3
import os
import shutil
import logging
import numpy as np
from osgeo import gdal
from datetime import datetime
import tempfile
from pystac_client import Client
import numexpr as ne

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(name)s - %(message)s")
logger = logging.getLogger("boreal_ndvi_trend")
# maap = MAAP()
# GDAL configurations used to successfully access LP DAAC Cloud Assets via vsicurl 
gdal.SetConfigOption('GDAL_HTTP_COOKIEFILE','~/cookies.txt')
gdal.SetConfigOption('GDAL_HTTP_COOKIEJAR', '~/cookies.txt')
gdal.SetConfigOption('GDAL_DISABLE_READDIR_ON_OPEN','EMPTY_DIR')
gdal.SetConfigOption('CPL_VSIL_CURL_ALLOWED_EXTENSIONS','TIF')
gdal.SetConfigOption('GDAL_HTTP_UNSAFESSL', 'YES')
gdal.SetConfigOption('GDAL_HTTP_MAX_RETRY', '10')
gdal.SetConfigOption('GDAL_HTTP_RETRY_DELAY', '0.5')
gdal.SetConfigOption('GDAL_HTTP_RETRY_DELAY', '0.5')
gdal.SetConfigOption('AWS_REQUEST_PAYER', 'requester')
os.environ["AWS_REQUEST_PAYER"] = "requester"

def split_s3_path(s3_path):
	path_parts=s3_path.replace("s3://","").split("/")
	bucket=path_parts.pop(0)
	key="/".join(path_parts)
	return bucket, key
def s3_key_exists(client, bucket, key):
	response = client.list_objects_v2(Bucket=bucket,Prefix=key)
	return 'Contents' in response
def tile_bounds(tile):
	xmin,x_sign,ymax,y_sign = int(tile[0:3]),tile[3],int(tile[4:6]),tile[6]
	if x_sign == 'W':
		xmin = -1*xmin
	if y_sign == 'S':
		ymax = -1*ymax
	xmax = xmin+1.0
	ymin = ymax-1.0
	return (xmin,ymin,xmax,ymax)
def mask_qa(qa_arr, mask_list=['dilated_cloud','cloud','shadow','snow','water']):
	QA_BIT = {'fill': 0,'dilated_cloud': 1,'cloud': 3,'shadow':4,'snow':5,'water':7}
	msk = np.zeros_like(qa_arr)
	for m in mask_list:
		if m in QA_BIT.keys():
			msk += (qa_arr & (1 << QA_BIT[m])) > 0
	return msk > 0
def index_setup(tile,indices,output):
	indices = indices.split(',')
	client = boto3.client('s3')
	index_bands = {'ndvi':['nir08','red'],'evi':['nir08','red','blue'],'savi':['nir08','red'],'msavi':['nir08','red'],'ndwi':['green','nir08'],'ndmi':['nir08','swir16'],'nbr2':['swir16','swir22'],'tvi':['nir08','red','green']}
	index_equations = {'ndvi':'(nir08-red)/(nir08+red)','evi':'2.5*(nir08-red)/(nir08+6*red-7.5*blue+1)','savi':'1.5*(nir08-red)/(nir08+red+0.5)','msavi':'0.5*(2*nir08+1-sqrt((2*nir08+1)*(2*nir08+1)))','ndwi':'(green-nir08)/(green+nir08)','ndmi':'(nir08-swir16)/(nir08+swir16)','nbr2':'(swir16-swir22)/(swir16+swir22)','tvi':'0.5*(120*(nir08-green)-200*(red-green))'}
	xmin,ymin,xmax,ymax = tile_bounds(tile)
	si_max = 10000
	si_min = -10000
	scaler = 10000
	nodata =-19999
	warp_options = gdal.WarpOptions(xRes=0.00025,yRes=0.00025,dstNodata=0,dstSRS='EPSG:4326',outputBounds=[xmin,ymin,xmax,ymax],format='COG',targetAlignedPixels=True,creationOptions=['COMPRESS=LZW'])
	qa_warp_options = gdal.WarpOptions(xRes=0.00025,yRes=0.00025,dstSRS='EPSG:4326',outputBounds=[xmin,ymin,xmax,ymax],format='COG',targetAlignedPixels=True,creationOptions=['COMPRESS=LZW'])
	url = 'https://landsatlook.usgs.gov/stac-server'
	collections = ['landsat-c2l2-sr']
	equations = {index: index_equations[index] for index in indices if index in index_equations}
	size = 4000
	sf = 0.0000275
	offset = -0.2
	val_min = 7273
	val_max = 43636
	bands = list(set([b for index in indices for b in index_bands[index]])) 
	
	for year in range(1984,2013):
		with tempfile.TemporaryDirectory() as tmpdir:
			catalog = Client.open(url)
			query = catalog.search(collections=collections, datetime=f'{year}-07/{year}-08',limit=50, bbox=(xmin-0.0125,ymin-0.0125,xmax+0.0125,ymax+0.015))
			for page_dict in query.pages_as_dicts():
				for item_dict in page_dict['features']:
					id = item_dict['id']
					input_array = np.full((len(bands),size,size),np.nan,dtype=np.float32)
					gt = None
					wkt = None
					mask_array = None
					process_flag = True
					date_str = ''
					try:
						date_str = id.split('_')[3]
						qa_file = item_dict['assets']['qa_pixel']['alternate']['s3']['href'].replace('s3://','/vsis3/').replace('usgs-landsat-ard/','usgs-landsat/')
						qa_warped = gdal.Warp(f'{tmpdir}/qa.tif', qa_file, options=qa_warp_options)
						qa_warped = None
						qa_ds = gdal.Open(f'{tmpdir}/qa.tif')
						gt = qa_ds.GetGeoTransform()
						wkt = qa_ds.GetProjectionRef()
						qa_band = qa_ds.GetRasterBand(1)
						qa_array = np.array(qa_band.ReadAsArray())
						mask_array = mask_qa(qa_array)
						qa_ds = None
						for b,band in enumerate(bands):
							if band in item_dict['assets']:
								band_file = item_dict['assets'][band]['alternate']['s3']['href']
							vsi_file = band_file.replace('s3://','/vsis3/')
							band_bucket,band_key=split_s3_path(band_file)
							if 'usgs-landsat-ard' not in band_file:
								band_warped = gdal.Warp(f'{tmpdir}/band.tif', vsi_file, options=warp_options)
								band_warped = None
								ds = gdal.Open(f'{tmpdir}/band.tif')
								gt = ds.GetGeoTransform()
								wkt = ds.GetProjectionRef()
								xsize, ysize=ds.RasterXSize,ds.RasterYSize
								band = ds.GetRasterBand(1)
								band_array = np.array(band.ReadAsArray())
								ds = None
								band_array_float = np.where((band_array>=val_min)&(band_array<=val_max),sf*band_array+offset,np.nan)
								band_array = None
								input_array[b,:,:] = band_array_float
								band_array_float = None
							else:
								logger.info(f'{band_file} not found')
								process_flag = False
					except Exception as err:
						logger.info(f'Loading QA error: {err}')
						process_flag = False
					if process_flag:
						input_arrays = {band:input_array[b,:,:] for b,band in enumerate(bands)}
						for si in equations:
							output_array = ne.evaluate(equations[si],local_dict=input_arrays)
							output_array = np.where(np.isnan(output_array),nodata,scaler*output_array)
							output_array = np.where(mask_array==1,nodata,output_array)
							output_array = np.where(output_array>si_max,nodata,output_array)
							output_array = np.where(output_array<si_min,nodata,output_array)
							output_array = output_array.astype(np.int16)
							driver = gdal.GetDriverByName('GTiff')
							output_ds = driver.Create(f'{tmpdir}/si.tif', size, size, 1, gdal.GDT_Int16, options=['COMPRESS=LZW','TILED=YES','COPY_SRC_OVERVIEWS=YES'])
							output_ds.SetGeoTransform(gt)
							output_ds.SetProjection(wkt)
							band = output_ds.GetRasterBand(1)
							band.WriteArray(output_array)
							band.SetNoDataValue(nodata)
							band.FlushCache()
							del output_ds
							output_bucket,output_key = split_s3_path(output)
							shutil.copyfile(f'{tmpdir}/si.tif',f'{output}/{tile}_{year}_{id}_{si}.tif')

		
def main():
	import argparse
	parser = argparse.ArgumentParser()
	parser.add_argument('--tile',type=str,required=True)
	parser.add_argument('--indices',type=str,required=False,default='ndvi,evi,savi,msavi,ndwi,ndmi,nbr2,tvi')
	parser.add_argument('--output',type=str,required=True)
	args = parser.parse_args()
	index_setup(args.tile,args.indices,args.output)
if __name__ == '__main__':
	main()
