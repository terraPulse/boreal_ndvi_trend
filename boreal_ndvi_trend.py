import numpy as np
from osgeo import gdal
import pandas as pd
from datetime import datetime
import time
from scipy import stats
import numpy as np
import earthaccess
import tempfile
import logging
import shutil
import os
os.environ['AWS_REQUEST_PAYER'] = 'requester'
os.environ['CPL_VSIL_CURL_NON_CACHED'] = '/vsis3/'
os.environ['GDAL_MAX_DATASET_POOL_SIZE'] = '4000'
# from maap.maap import MAAP
# from boto3 import Session
import boto3
import requests

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(name)s - %(message)s")
logger = logging.getLogger("boreal_si_trend")
# maap = MAAP()
# GDAL configurations used to successfully access LP DAAC Cloud Assets via vsicurl 
gdal.SetConfigOption('GDAL_HTTP_COOKIEFILE','~/cookies.txt')
gdal.SetConfigOption('GDAL_HTTP_COOKIEJAR', '~/cookies.txt')
gdal.SetConfigOption('GDAL_DISABLE_READDIR_ON_OPEN','EMPTY_DIR')
gdal.SetConfigOption('CPL_VSIL_CURL_ALLOWED_EXTENSIONS','TIF')
gdal.SetConfigOption('GDAL_HTTP_UNSAFESSL', 'YES')
gdal.SetConfigOption('GDAL_HTTP_MAX_RETRY', '10')
gdal.SetConfigOption('GDAL_HTTP_RETRY_DELAY', '0.5')
gdal.SetConfigOption('AWS_REQUEST_PAYER', 'requester')

def split_s3_path(s3_path):
    path_parts=s3_path.replace("s3://","").split("/")
    bucket=path_parts.pop(0)
    key="/".join(path_parts)
    return bucket, key
def s3_key_exists(client, bucket, key):
    response = client.list_objects_v2(Bucket=bucket,Prefix=key)
    return 'Contents' in response
def mask_hls(qa_arr, mask_list=['cloud','adj_cloud','water','snowice']):
	QA_BIT = {'cirrus': 0,'cloud': 1,'adj_cloud': 2,'cloud shadow':3,'snowice':4,'water':5,'aerosol_l': 6,'aerosol_h': 7}
	msk = np.zeros_like(qa_arr)
	for m in mask_list:
		if m in QA_BIT.keys():
			msk += (qa_arr & (1 << QA_BIT[m])) > 0
		if m == 'aerosol_high':
			msk += ((qa_arr & (1 << QA_BIT['aerosol_h'])) > 0) * ((qa_arr & (1 << QA_BIT['aerosol_l'])) > 0)
		if m == 'aerosol_moderate':
			msk += ((qa_arr & (1 << QA_BIT['aerosol_h'])) > 0) * ((qa_arr | (1 << QA_BIT['aerosol_l'])) != qa_arr)
		if m == 'aerosol_low':
			msk += ((qa_arr | (1 << QA_BIT['aerosol_h'])) != qa_arr) * ((qa_arr & (1 << QA_BIT['aerosol_l'])) > 0)
	return msk > 0

def tile_bounds(tile):
	xmin,x_sign,ymax,y_sign = int(tile[0:3]),tile[3],int(tile[4:6]),tile[6]
	if x_sign == 'W':
		xmin = -1*xmin
	if y_sign == 'S':
		ymax = -1*ymax
	xmax = xmin+1.0
	ymin = ymax-1.0
	return (xmin,ymin,xmax,ymax)
def get_files(bbox, start_date, end_date,si):
	file_list = []
	results = earthaccess.search_data(short_name='HLSL30_VI',cloud_hosted=True,temporal = (start_date, end_date),bounding_box = bbox)
	for rec in results:
		for url in rec.data_links(access='direct'):
			if f'{si.upper()}.tif' in url:
				file_list.append(url)
	results = earthaccess.search_data(short_name='HLSS30_VI',cloud_hosted=True,temporal = (start_date, end_date),bounding_box = bbox)
	for rec in results:
		for url in rec.data_links(access='direct'):
			if f'{si.upper()}.tif' in url:
				file_list.append(url)
	return file_list
def get_files_lc2(tile,year,si):
	s3_client = boto3.client('s3')
	response = s3_client.get_object(Bucket='maap-ops-workspace', Key=f'pswang/dps_output/boreal_si/{si}_list.csv')
	df = pd.read_csv(response['Body'])
	df_tile = df[(df['tile'] == tile) & (df['year'] == year)]
	prod_ids = list(set(df_tile['prod_id'].tolist()))
	file_list = []
	for prod_id in prod_ids:
		file_list.append(df_tile[(df_tile['prod_id'] == prod_id)]['file'].tolist()[0])
	return file_list

def boreal_si_trend(tile,ys,ye,ds,de,output,indices,mode):
	bbox = tile_bounds(tile)
	nodata = -19999
	warp_options = gdal.WarpOptions(xRes=0.00025,yRes=0.00025,dstSRS='EPSG:4326',outputBounds=bbox,targetAlignedPixels=True)
	indices = indices.split(',')
	for si in indices:
		gt = None
		wkt = None
		xsize = None
		ysize = None
		n = None
		sum_x = None
		sum_y = None
		sum_xx = None
		sum_yy = None
		sum_xy = None
		with tempfile.TemporaryDirectory() as tmpdir:
			for year in range(ys,ye+1):
				if mode == 'hls':
					si_files = get_files(bbox,f'{year}-{ds}',f'{year}-{de}',si)
				elif mode == 'lc2':
					si_files = get_files_lc2(tile,year,si)
				with tempfile.TemporaryDirectory() as cachedir:
					for si_file in si_files:
						logger.info(si_file)
						si_arr = None
						try:
							if mode == 'hls':
								local_file = earthaccess.download(si_file,local_path=cachedir,provider='LPCLOUD')
								logger.info(local_file)
								si_warped = gdal.Warp(f'{tmpdir}/si.tif', local_file, options=warp_options)
								si_warped = None
								si_ds = gdal.Open(f'{tmpdir}/si.tif')
								si_arr = si_ds.GetRasterBand(1).ReadAsArray()
								gt = si_ds.GetGeoTransform()
								wkt = si_ds.GetProjectionRef()
								xsize, ysize=si_ds.RasterXSize,si_ds.RasterYSize
								si_ds = None
								local_file_fmask = earthaccess.download(si_file.replace('_VI','').replace('-VI','').replace(f'{si.upper()}.tif','Fmask.tif'),local_path=cachedir,provider='LPCLOUD')
								logger.info(local_file_fmask)
								qa_warped = gdal.Warp(f'{tmpdir}/qa.tif', local_file_fmask, options=warp_options)
								qa_warped = None
								qa_ds = gdal.Open(f'{tmpdir}/qa.tif')
								qa_arr = qa_ds.GetRasterBand(1).ReadAsArray()
								qa_ds = None
								mask = mask_hls(qa_arr,mask_list=['cloud','adj_cloud','water','snowice'])
								si_arr = np.where((si_arr == nodata) |(mask == 1),np.nan,0.0001*si_arr)
							elif mode == 'lc2':
								si_ds = gdal.Open(si_file)
								si_arr = si_ds.GetRasterBand(1).ReadAsArray()
								gt = si_ds.GetGeoTransform()
								wkt = si_ds.GetProjectionRef()
								xsize, ysize=si_ds.RasterXSize,si_ds.RasterYSize
								si_ds = None
								si_arr = np.where(si_arr == nodata,np.nan,0.0001*si_arr)
							if n is None:
								n = np.where(np.isnan(si_arr),0,1)
								sum_x = np.where(np.isnan(si_arr),0,year)
								sum_y = np.where(np.isnan(si_arr),0,si_arr)
								sum_xx = np.where(np.isnan(si_arr),0,year*year)
								sum_yy = np.where(np.isnan(si_arr),0,si_arr*si_arr)
								sum_xy = np.where(np.isnan(si_arr),0,year*si_arr)
							else:
								n += np.where(np.isnan(si_arr),0,1)
								sum_x += np.where(np.isnan(si_arr),0,year)
								sum_y += np.where(np.isnan(si_arr),0,si_arr)
								sum_xx += np.where(np.isnan(si_arr),0,year*year)
								sum_yy += np.where(np.isnan(si_arr),0,si_arr*si_arr)
								sum_xy += np.where(np.isnan(si_arr),0,year*si_arr)
						except Exception as e:
							logger.info(f'Error processing {si_file}: {e}')
			with np.errstate(divide='ignore', invalid='ignore'):
				denom = (n * sum_xx) - (sum_x**2)
				slopes = ((n * sum_xy) - (sum_x * sum_y)) / denom
				intercepts = (sum_y - slopes * sum_x) / n
				s_xx = sum_xx - (sum_x**2 / n)
				s_yy = sum_yy - (sum_y**2 / n)
				s_xy = sum_xy - (sum_x * sum_y / n)
				rss = np.maximum(s_yy - (s_xy**2 / s_xx), 0)
				df = n - 2
				s_err = np.sqrt(rss / df)
				se_slope = s_err / np.sqrt(s_xx)
				t_stats = slopes / se_slope
				p_values = stats.t.sf(np.abs(t_stats), df) * 2
				slopes[(n < 2) | np.isnan(slopes)] = nodata
				intercepts[(n < 2) | np.isnan(intercepts)] = nodata
				p_values[(n < 2) | np.isnan(p_values)] = nodata
	
				driver = gdal.GetDriverByName('GTiff')
				output_ds = driver.Create(f'{tmpdir}/{tile}_slope.tif', xsize,ysize, 1, gdal.GDT_Float32, options=['COMPRESS=LZW','TILED=YES','COPY_SRC_OVERVIEWS=YES','BIGTIFF=YES'])
				output_ds.SetGeoTransform(gt)
				output_ds.SetProjection(wkt)
				band = output_ds.GetRasterBand(1)
				band.WriteArray(slopes)
				band.SetNoDataValue(nodata)
				band.FlushCache()
				del output_ds
				driver = gdal.GetDriverByName('GTiff')
				output_ds = driver.Create(f'{tmpdir}/{tile}_intercept.tif', xsize,ysize, 1, gdal.GDT_Float32, options=['COMPRESS=LZW','TILED=YES','COPY_SRC_OVERVIEWS=YES','BIGTIFF=YES'])
				output_ds.SetGeoTransform(gt)
				output_ds.SetProjection(wkt)
				band = output_ds.GetRasterBand(1)
				band.WriteArray(intercepts)
				band.SetNoDataValue(nodata)
				band.FlushCache()
				del output_ds
				driver = gdal.GetDriverByName('GTiff')
				output_ds = driver.Create(f'{tmpdir}/{tile}_pval.tif', xsize,ysize, 1, gdal.GDT_Float32, options=['COMPRESS=LZW','TILED=YES','COPY_SRC_OVERVIEWS=YES','BIGTIFF=YES'])
				output_ds.SetGeoTransform(gt)
				output_ds.SetProjection(wkt)
				band = output_ds.GetRasterBand(1)
				band.WriteArray(p_values)
				band.SetNoDataValue(nodata)
				band.FlushCache()
				del output_ds
				n = n.astype('uint16')
				driver = gdal.GetDriverByName('GTiff')
				output_ds = driver.Create(f'{tmpdir}/{tile}_num.tif', xsize,ysize, 1, gdal.GDT_UInt16, options=['COMPRESS=LZW','TILED=YES','COPY_SRC_OVERVIEWS=YES','BIGTIFF=YES'])
				output_ds.SetGeoTransform(gt)
				output_ds.SetProjection(wkt)
				band = output_ds.GetRasterBand(1)
				band.WriteArray(n)
				band.SetNoDataValue(0)
				band.FlushCache()
				del output_ds
				logger.info("Copying files")
				shutil.copyfile(f'{tmpdir}/{tile}_slope.tif',f'{output}/{tile}_{mode}_{si}_slope.tif')
				shutil.copyfile(f'{tmpdir}/{tile}_intercept.tif',f'{output}/{tile}_{mode}_{si}_intercept.tif')
				shutil.copyfile(f'{tmpdir}/{tile}_pval.tif',f'{output}/{tile}_{mode}_{si}_pval.tif')
				shutil.copyfile(f'{tmpdir}/{tile}_num.tif',f'{output}/{tile}_{mode}_{si}_num.tif')
			

def main():
	import argparse
	parser = argparse.ArgumentParser()
	parser.add_argument('--tile',type=str,required=True)
	parser.add_argument('--ys',type=int,required=True)
	parser.add_argument('--ye',type=int,required=True)
	parser.add_argument('--ds',type=str,required=True)
	parser.add_argument('--de',type=str,required=True)
	parser.add_argument('--output',type=str,required=True)
	parser.add_argument('--indices',type=str,required=False,default='ndvi,evi,savi,msavi,ndwi,ndmi,nbr2,tvi')
	parser.add_argument('--mode',type=str,required=True)
	parser.add_argument('--user',type=str,required=True)
	parser.add_argument('--pwd',type=str,required=True)
	# parser.add_argument('--token',type=str,required=True)
	
	args = parser.parse_args()
	os.environ["EARTHDATA_USERNAME"] = args.user
	os.environ["EARTHDATA_PASSWORD"] = args.pwd
	logger.info("Username: "+os.environ.get('EARTHDATA_USERNAME'))
	logger.info("Password: "+os.environ.get('EARTHDATA_PASSWORD'))
	
	earthaccess.login()
	boreal_si_trend(args.tile,args.ys,args.ye,args.ds,args.de,args.output,args.indices,args.mode)
if __name__ == '__main__':
	main()
