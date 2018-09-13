
//#define ITEM_COUNT 2
#define _PI 3.14159265358979323846
#define _PI2 1.57079632679489661923
#define _RAD 6372795



#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <math.h>
#include <iostream>
#include <stdio.h>
#include <string>

using namespace std;

cudaError_t addWithCuda(int *c, const int *a, const int *b, unsigned int size);


__global__ void geo_invert(double2* d_dot1, double2* d_dot2, double* d_dist, double* d_azimut, long count)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	if (idx < count)
	{
		d_dot1[idx].x = d_dot1[idx].x * _PI / 180;	//lat1
		d_dot1[idx].y = d_dot1[idx].y * _PI / 180;	//lng1
		d_dot2[idx].x = d_dot2[idx].x * _PI / 180;	//lat2
		d_dot2[idx].y = d_dot2[idx].y * _PI / 180;	//lng2

		double cl1, cl2, sl1, sl2, delta, cdelta, sdelta;
		cl1 = cos(d_dot1[idx].x);
		cl2 = cos(d_dot2[idx].x);
		sl1 = sin(d_dot1[idx].x);
		sl2 = sin(d_dot2[idx].x);
		delta = d_dot2[idx].y - d_dot1[idx].y;
		cdelta = cos(delta);
		sdelta = sin(delta);

		double x, y, z, ad, z2;
		y = sqrt(pow(cl2*sdelta, 2) + pow(cl1*sl2 - sl1*cl2*cdelta, 2));
		x = sl1*sl2 + cl1*cl2*cdelta;
		ad = atan(y / x);
		d_dist[idx] = ad * _RAD;

		x = (cl1*sl2) - (sl1*cl2*cdelta);
		y = sdelta*cl2;

		if (x == 0)
		{
			if (y > 0)
				z = -90;
			else if (y < 0)
				z = 90;
			else if (y == 0)
				z = 0;
		}
		else
		{
			z = atan(-y / x) * 180 / _PI;
			if (x < 0)
			{
				z = z + 180;
			}
		}

		z2 = z + 180.0f;

		while (z2 >= 360)
		{
			z2 = z2 - 360;
		}

		z2 = z2 - 180;


		z2 = -z2 * _PI / 180;
		double anglerad2;
		anglerad2 = z2 - ((2 * _PI) * floor(z2 / (2 * _PI)));
		d_azimut[idx] = anglerad2 * 180 / _PI;


	}
}

__device__ double CartToSpher(double3* x, double2* y)
{
	double p;	

	p = hypot(x->x, x->y); //0,566796731779912 (0,566796731779913) 0,000000000000001
	y->y = atan(x->y / x->x); //0,658744870833875 (0,658744870833875) 0
	y->x = atan(x->z / p); //0,968183828701654 (0,968183828701654) 0

	return hypot(p, x->z);
}

__device__ void SpherToCart(double2* y, double3* x)
{
	double p;

	p = cos(y->x); //0,509535037739044 (0,509535037739044) 0
	x->z = sin(y->x); //0,860449908661899 (0,860449908661899)  0
	x->y = p * sin(y->y); //sin -0,896141078848377 () ; -0,456615278430515379561368931588 (-0,456615278430516) -0,0000000000000007
	x->x = p * cos(y->y); //cos -0,443769272032738 () ; -0,226115992772629239921468822472 (-0,226115992772629)  0,0000000000000002

	return;
}

__device__ void Rotate2(double3* x, double a)
{
	double c, s, xj;

	c = cos(a); //-0,0246449569641315 (-0,0246449569641315) 0
	s = sin(a); //-0,999696266921226 (-0,999696266921226) 0
	xj = x->x * c + x->y * s; //0,448200835371488472286581966712384075135922396544509757722076 (0,448200835371489) 0,0000000000000006
	x->y = -x->x * s + x->y * c; //0,346950351388623496054477853457869740953054474470049783262896 (0,346950351388624) 0,0000000000000006
	x->x = xj; //0,448200835371488472286581966712384075135922396544509757722076 (0,448200835371489) 0,0000000000000006

	return;
}

__device__ void Rotate1(double3* x, double a)
{
	double c, s, xj;

	c = cos(a); //0,799692643650457 (0,799692643650457) 0
	s = sin(a); //-0,60040958993952 (-0,60040958993952) 0
	xj = x->z * c + x->x * s; //0,82385767268601006086792194931074587209689344 (0,823857672686011) 0,000000000000001
	x->x = -x->z * s + x->x * c; //0,335799080791196711492751735883669680645330296 (0,335799080791197) 0,0000000000000003
	x->z = xj; //0,82385767268601006086792194931074587209689344 (0,823857672686011) 0,000000000000001

	return;
}

__device__ void SphereDirect(double2* pt1, double azi, double dist, double2* pt2)
{
	double2 pt;
	double3 x;

	pt.x = _PI2 - dist; //1,036151994127676035482994172886465044 (1,03615199412768)	0,000000000000004 15
	pt.y = _PI - azi;  // -2,030596755794815247590950216888888888888888888888888889 (-2,03059675579482)  -0,000000000000005 15

	SpherToCart(&pt, &x);               // сферические -> декартовы
	Rotate1(&x, pt1->x - _PI2); // первое вращение
	Rotate2(&x, -pt1->y);           // второе вращение
	CartToSpher(&x, pt2);           // декартовы -> сферические 

}

__global__ void geo_direct(double2* d_dot1, double* d_dist, double* d_azimut, double2* d_dot2, long count)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	if (idx < count)
	{
		d_dot1[idx].x = d_dot1[idx].x * _PI / 180;	//pt1[0] 0,926783132217889 (0.926783132217889) 0
		d_dot1[idx].y = d_dot1[idx].y * _PI / 180;	//pt1[1] 1,595443779225097 (1,5954437792251) 0,000000000000003 15
													//d_azimut[idx] 296,344624 (296.344624)
		d_azimut[idx] = d_azimut[idx] * _PI / 180; //5,172189409384608486050950216 (5,17218940938461) 0,000000000000002 15
		d_dist[idx] = d_dist[idx] / _RAD; //0,534644332667220583747005827113534956 (0,534644332667221)	0,0000000000000005 15

		double2 pt2;

		SphereDirect(&d_dot1[idx], d_azimut[idx], d_dist[idx], &pt2);

		d_dot2[idx].x = pt2.x * 180 / _PI; //pt 0,968183828701654 (0,968183734421639) // 55,472847177421830384141692371164405517504840227 (55,4728417755749) // -0,000000094280015 // -0,00000540184693
		d_dot2[idx].y = pt2.y * 180 / _PI; //pt 0,658744870833875 (0,658744865257753) // 37,743300874671594906787964611067062014639516391 (37,7433005551833) // -0,000000005576122 // -0,00000031948829

		if (d_dot2[idx].x < 0)
			d_dot2[idx].x += 180;

		if (d_dot2[idx].y < 0)
			d_dot2[idx].y += 180;

	}
}

/*__global__ void geo_invert(double2* d_dot1, double2* d_dot2, double* d_dist, double* d_azimut, long count)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	if (idx < count)
	{
		double3 x;
		double2 pt;

		d_dot1[idx].x = d_dot1[idx].x * _PI / 180;	//lat1
		d_dot1[idx].y = d_dot1[idx].y * _PI / 180;	//lng1
		d_dot2[idx].x = d_dot2[idx].x * _PI / 180;	//lat2
		d_dot2[idx].y = d_dot2[idx].y * _PI / 180;	//lng2


		SpherToCart(&d_dot2[idx], &x);			// сферические -> декартовы
		Rotate2(&x, d_dot1[idx].y);			// первое вращение
		Rotate1(&x, _PI2 - d_dot1[idx].x);	// второе вращение
		CartToSpher(&x, &pt);	     		// декартовы -> сферические
		d_azimut[idx] = (_PI - pt.y)* 180 / _PI;
		d_dist[idx] = (_PI2 - pt.x) * _RAD;
	}
}*/

__device__ double d_abs(double var)
{
	if (var<0)
	{
		return -var;
	}
	else
	{
		return var;
	}
}

__global__ void d_cudainit(int *a, int *b)
{
    int i = threadIdx.x;
	if (i==1)
	{
		b[i] = a[i] * 2;
	}
}

__device__ void SpherToCartR(double* y, double* x)
{
	double y_lat = y[0] * _PI / 180;
	double y_lon = y[1] * _PI / 180;

	double cos_y_lat = cos(y_lat);

	x[0] = _RAD * cos_y_lat * cos(y_lon);
	x[1] = _RAD * cos_y_lat * sin(y_lon);
	x[2] = _RAD * sin(y_lat);

	return;
}

__device__ bool LineIsTooFar(double* M, double* A, double* B, float max_delta)
{
	double lat_lag = 0.00001 * max_delta * 1.2;
	double lon_lag = 0.00003 * max_delta * 1.2;

	double sqr_left_low_lat = M[0] - lat_lag;
	double sqr_left_low_lon = M[1] - lon_lag;
	double sqr_right_high_lat = M[0] + lat_lag;
	double sqr_right_high_lon = M[1] + lon_lag;

	if (A[1] < sqr_left_low_lon && B[1] < sqr_left_low_lon)
		return true;

	if (A[1] > sqr_right_high_lon && B[1] > sqr_right_high_lon)
		return true;

	if (A[0] < sqr_left_low_lat && B[0] < sqr_left_low_lat)
		return true;

	if (A[0] > sqr_right_high_lat && B[0] > sqr_right_high_lat)
		return true;

	return false;
}

__device__ double LineLength(double* a, double* b)
{
	double K1, K2, K3, K4;

	K1 = b[0] * a[0] + b[1] * a[1] + b[2] * a[2];
	K2 = sqrt(pow(b[0], 2) + pow(b[1], 2) + pow(b[2], 2));
	K3 = sqrt(pow(a[0], 2) + pow(a[1], 2) + pow(a[2], 2));
	K4 = K1 / (K2 * K3);

	if (K4 > 1) K4 = 1;

	return _RAD * acos(K4);
}

__device__ bool M_ProjectionOnPlane(double* m_dec, double plane_A, double plane_B, double plane_C, double* a_dec, double* b_dec)
{
	double t = -(plane_A * m_dec[0] + plane_B * m_dec[1] + plane_C * m_dec[2]) / (pow(plane_A, 2) + pow(plane_B, 2) + pow(plane_C, 2));

	double* dot_on_plane = new double[3];
	dot_on_plane[0] = plane_A * t + m_dec[0];
	dot_on_plane[1] = plane_B * t + m_dec[1];
	dot_on_plane[2] = plane_C * t + m_dec[2];

	double* dot_k = new double[3];
	double K = sqrt(pow(dot_on_plane[0], 2) + pow(dot_on_plane[1], 2) + pow(dot_on_plane[2], 2));
	dot_k[0] = (_RAD * dot_on_plane[0]) / K;
	dot_k[1] = (_RAD * dot_on_plane[1]) / K;
	dot_k[2] = (_RAD * dot_on_plane[2]) / K;

	double line_AB = LineLength(a_dec, b_dec);
	double line_AK = LineLength(a_dec, dot_k);
	double line_BK = LineLength(b_dec, dot_k);

	delete dot_on_plane;
	delete dot_k;

	return d_abs(line_AK + line_BK - line_AB) < 0.01;
}

__device__ double DistanceToLine(double* m_dec, double* a_dec, double* b_dec, bool calc_mb)
{
	if (a_dec[0] == b_dec[0] &&
		a_dec[1] == b_dec[1] &&
		a_dec[2] == b_dec[2])
	{
		return LineLength(m_dec, a_dec);
	}

	double plane_A, plane_B, plane_C;
	plane_A = a_dec[1] * b_dec[2] - a_dec[2] * b_dec[1];
	plane_B = a_dec[2] * b_dec[0] - a_dec[0] * b_dec[2];
	plane_C = a_dec[0] * b_dec[1] - a_dec[1] * b_dec[0];

	double d, MK_length, MA_length, MB_length, minlength;

	if (M_ProjectionOnPlane(m_dec, plane_A, plane_B, plane_C, a_dec, b_dec))
	{
		d = d_abs(plane_A*m_dec[0] + plane_B * m_dec[1] + plane_C*m_dec[2]) / sqrt(pow(plane_A, 2) + pow(plane_B, 2) + pow(plane_C, 2));
		MK_length = _RAD * asin(d / _RAD);
		minlength = MK_length;
	}
	else
	{
		MA_length = LineLength(m_dec, a_dec);
		if (calc_mb)
		{
			MB_length = LineLength(m_dec, b_dec);
			minlength = fmin(MA_length, MB_length);
		}
		else
		{
			minlength = MA_length;
		}
	}

	return minlength;
}

__device__ int dot_near_polyline(double dot_lat, double dot_lon, double* line_lat, double* line_lon, long line_count, float max_delta)
{
	double* m_sph = new double[2];
	m_sph[0] = dot_lat;
	m_sph[1] = dot_lon;

	double* m_dec = new double[3];

	SpherToCartR(m_sph, m_dec);

	double disttoline;
	int dotisclose = 0;

	double* a_sph = new double[2];
	double* b_sph = new double[2];
	double* a_dec = new double[3];
	double* b_dec = new double[3];

	for (long i = 0; i <= line_count - 2; i++)
	{
		//line dots
		a_sph[0] = line_lat[i];
		a_sph[1] = line_lon[i];

		b_sph[0] = line_lat[i+1];
		b_sph[1] = line_lon[i+1];

		if (LineIsTooFar(m_sph, a_sph, b_sph, max_delta))
		{
			continue;
		}

		SpherToCartR(a_sph, a_dec);
		SpherToCartR(b_sph, b_dec);

		disttoline = DistanceToLine(m_dec, a_dec, b_dec, i == line_count - 2);

		if (disttoline < max_delta)
		{
			dotisclose = 1;
			break;
		}
	}

	delete m_sph;
	delete m_dec;
	delete a_sph;
	delete b_sph;
	delete a_dec;
	delete b_dec;

	return dotisclose;



}

__global__ void dotarray_near_polyline(double* dot_lat, double* dot_lon, double* line_lat, double* line_lon, long dot_count, long line_count, float max_delta, int* dot_result)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	if (idx < dot_count)
	{
		dot_result[idx] = dot_near_polyline(dot_lat[idx], dot_lon[idx], line_lat, line_lon, line_count, max_delta);
		//one thread calculate 1 dot near line array
	}
}


extern "C" __declspec(dllexport)	int GetInvertGeo(double* dot1_lat, double* dot1_lon, double* dot2_lat, double* dot2_lon, double* dist, double* azimut, long count)
{
	const int blockSize = 1024;
	int numOfBlocks = (count + blockSize - 1) / blockSize;
	dim3 dimGrid(numOfBlocks);
	dim3 dimBlock(blockSize);
	cudaError_t cudaStatus;

	int size_double = sizeof(double) * count;
	int size_double2 = sizeof(double2) * count;

	double2 *d_dot1, *d_dot2;
	double *d_azimut;
	double *d_dist;
	cudaStatus = cudaMalloc((void**)&d_dot1, size_double2);
	if (cudaStatus != cudaSuccess)
	{
		return 2;
	}
	cudaStatus = cudaMalloc((void**)&d_dot2, size_double2);
	if (cudaStatus != cudaSuccess)
	{
		return 2;
	}
	cudaStatus = cudaMalloc((void**)&d_azimut, size_double);
	if (cudaStatus != cudaSuccess)
	{
		return 2;
	}
	cudaStatus = cudaMalloc((void**)&d_dist, size_double);
	if (cudaStatus != cudaSuccess)
	{
		return 2;
	}

	double2 *dot1 = new double2[count];
	double2 *dot2 = new double2[count];

	for (long i = 0; i < count; i++)
	{
		dot1[i].x = dot1_lat[i];
		dot1[i].y = dot1_lon[i];
		dot2[i].x = dot2_lat[i];
		dot2[i].y = dot2_lon[i];
	}


	//перенести входные массивы в видеопамять
	cudaStatus = cudaMemcpy(d_dot1, dot1, size_double2, cudaMemcpyKind::cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess)
	{
		return 3;
	}
	cudaStatus = cudaMemcpy(d_dot2, dot2, size_double2, cudaMemcpyKind::cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess)
	{
		return 3;
	}

	delete dot1;
	delete dot2;

	geo_invert <<< dimGrid, dimBlock >>> (d_dot1, d_dot2, d_dist, d_azimut, count);
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		return 4;
	}

	//double *dist = new double[ITEM_COUNT];
	//double *azimut = new double[ITEM_COUNT];
	cudaStatus = cudaMemcpy(azimut, d_azimut, size_double, cudaMemcpyKind::cudaMemcpyDeviceToHost);
	if (cudaStatus != cudaSuccess)
	{
		return 3;
	}
	cudaStatus = cudaMemcpy(dist, d_dist, size_double, cudaMemcpyKind::cudaMemcpyDeviceToHost);
	if (cudaStatus != cudaSuccess)
	{
		return 3;
	}


	cudaStatus = cudaFree(d_dot1);
	if (cudaStatus != cudaSuccess)
	{
		return 5;
	}
	cudaStatus = cudaFree(d_dot2);
	if (cudaStatus != cudaSuccess)
	{
		return 5;
	}
	cudaStatus = cudaFree(d_azimut);
	if (cudaStatus != cudaSuccess)
	{
		return 5;
	}
	cudaStatus = cudaFree(d_dist);
	if (cudaStatus != cudaSuccess)
	{
		return 5;
	}

	//вывод результатов
	//cout.precision(6);
	//cout << dist[0] << endl << azimut[0] << endl;

	//delete dist;
	//delete azimut;
	return 0;
}

extern "C" __declspec(dllexport)	int GetDirectGeo(double* dot1_lat, double* dot1_lon, double* dist, double* azimut,  double* dot2_lat, double* dot2_lon, long count)
{
	const int blockSize = 1024;
	//int numOfBlocks = (ITEM_COUNT + blockSize - 1) / blockSize;
	int numOfBlocks = (count - 1) / blockSize + 1;
	dim3 dimGrid(numOfBlocks);
	dim3 dimBlock(blockSize);
	cudaError_t cudaStatus;

	int size_double = sizeof(double) * count;
	int size_double2 = sizeof(double2) * count;

	double2 *d_dot1, *d_dot2;
	double *d_azimut;
	double *d_dist;
	cudaStatus = cudaMalloc((void**)&d_dot1, size_double2);
	if (cudaStatus != cudaSuccess)
	{
		return 2;
	}
	cudaStatus = cudaMalloc((void**)&d_dot2, size_double2);
	if (cudaStatus != cudaSuccess)
	{
		return 2;
	}
	cudaStatus = cudaMalloc((void**)&d_azimut, size_double);
	if (cudaStatus != cudaSuccess)
	{
		return 2;
	}
	cudaStatus = cudaMalloc((void**)&d_dist, size_double);
	if (cudaStatus != cudaSuccess)
	{
		return 2;
	}


	//массивы для хоста
	//double2 *dot1 = new double2[ITEM_COUNT];
	//double *azimut = new double[ITEM_COUNT];
	//double *dist = new double[ITEM_COUNT];

	//начальная точка
	//dot1[0].x = 53.100762;
	//dot1[0].y = 91.412195;
	//azimut[0] = 296.344624;
	//dist[0] = 3407178.73;


	double2 *dot1 = new double2[count];
	double2 *dot2 = new double2[count];

	for (long i = 0; i < count; i++)
	{
		dot1[i].x = dot1_lat[i];
		dot1[i].y = dot1_lon[i];
	}

	cudaStatus = cudaMemcpy(d_dot1, dot1, size_double2, cudaMemcpyKind::cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess)
	{
		return 3;
	}
	cudaStatus = cudaMemcpy(d_azimut, azimut, size_double, cudaMemcpyKind::cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess)
	{
		return 3;
	}
	cudaStatus = cudaMemcpy(d_dist, dist, size_double, cudaMemcpyKind::cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess)
	{
		return 3;
	}

	delete dot1;

	geo_direct <<< dimGrid, dimBlock >>> (d_dot1, d_dist, d_azimut,  d_dot2, count);
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		return 4;
	}


	//double2 *dot2 = new double2[ITEM_COUNT];
	cudaStatus = cudaMemcpy(dot2, d_dot2, size_double2, cudaMemcpyKind::cudaMemcpyDeviceToHost);
	if (cudaStatus != cudaSuccess)
	{
		return 3;
	}

	for (long i = 0; i < count; i++)
	{
		dot2_lat[i] = dot2[i].x;
		dot2_lon[i] = dot2[i].y;
	}

	delete dot2;

	cudaStatus = cudaFree(d_dot1);
	if (cudaStatus != cudaSuccess)
	{
		return 5;
	}
	cudaStatus = cudaFree(d_dot2);
	if (cudaStatus != cudaSuccess)
	{
		return 5;
	}
	cudaStatus = cudaFree(d_azimut);
	if (cudaStatus != cudaSuccess)
	{
		return 5;
	}
	cudaStatus = cudaFree(d_dist);
	if (cudaStatus != cudaSuccess)
	{
		return 5;
	}

	//вывод результатов
	//cout.precision(10);
	//cout << dot2[0].x << endl << dot2[0].y << endl;

	//delete dist;
	//delete azimut;
	return 0;
}

extern "C" __declspec(dllexport)	int CudaInitialize()
{
    int* a = new int[1];
    int* b = new int[1];
	cudaError_t cudaStatus;

	const int blockSize = 1024;
	int numOfBlocks = 1;
	dim3 dimGrid(numOfBlocks);
	dim3 dimBlock(blockSize);

	cudaStatus = cudaSetDevice(0);
	if (cudaStatus != cudaSuccess)
	{
		return 1;
	}

	int *d_a;
	int *d_b;
	cudaStatus = cudaMalloc((void**)&d_a, sizeof(int));
	if (cudaStatus != cudaSuccess)
	{
		return 2;
	}
	cudaStatus = cudaMalloc((void**)&d_b, sizeof(int));
	if (cudaStatus != cudaSuccess)
	{
		return 2;
	}

	a[0] = 1;
	
	cudaStatus = cudaMemcpy(d_a, a, sizeof(int), cudaMemcpyKind::cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess)
	{
		return 3;
	}

	d_cudainit <<< dimGrid, dimBlock >>> (d_a, d_b);
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		return 4;
	}

	cudaStatus = cudaMemcpy(b, d_b, sizeof(int), cudaMemcpyKind::cudaMemcpyDeviceToHost);
	if (cudaStatus != cudaSuccess)
	{
		return 3;
	}

	cudaStatus = cudaFree(d_a);
	if (cudaStatus != cudaSuccess)
	{
		return 5;
	}
	cudaStatus = cudaFree(d_b);
	if (cudaStatus != cudaSuccess)
	{
		return 5;
	}

	delete a;
	delete b;

    return 0;
}

extern "C" __declspec(dllexport)	int DotArrayNearPolyline(double* dot_lat, double* dot_lon, double* line_lat, double* line_lon, long dot_count, long line_count, float max_delta, int* dot_result, char str[])
{
	const int blockSize = 1024;
	//int numOfBlocks = (ITEM_COUNT + blockSize - 1) / blockSize;
	int numOfBlocks = (dot_count - 1) / blockSize + 1;
	dim3 dimGrid(numOfBlocks);
	dim3 dimBlock(blockSize);
	cudaError_t cudaStatus;

	int size_double_dots = sizeof(double) * dot_count;
	int size_double_polyline = sizeof(double) * line_count;
	int size_int = sizeof(int) * dot_count;

	double *d_dot_lat;
	double *d_dot_lon;
	double *d_line_lat;
	double *d_line_lon;
	int *d_dot_result;

	cudaStatus = cudaMalloc((void**)&d_dot_lat, size_double_dots);
	if (cudaStatus != cudaSuccess)
	{
		return 2;
	}
	cudaStatus = cudaMalloc((void**)&d_dot_lon, size_double_dots);
	if (cudaStatus != cudaSuccess)
	{
		return 2;
	}
	cudaStatus = cudaMalloc((void**)&d_line_lat, size_double_polyline);
	if (cudaStatus != cudaSuccess)
	{
		return 2;
	}
	cudaStatus = cudaMalloc((void**)&d_line_lon, size_double_polyline);
	if (cudaStatus != cudaSuccess)
	{
		return 2;
	}
	cudaStatus = cudaMalloc((void**)&d_dot_result, size_int);
	if (cudaStatus != cudaSuccess)
	{
		return 2;
	}

	cudaStatus = cudaMemcpy(d_dot_lat, dot_lat, size_double_dots, cudaMemcpyKind::cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess)
	{
		return 3;
	}
	cudaStatus = cudaMemcpy(d_dot_lon, dot_lon, size_double_dots, cudaMemcpyKind::cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess)
	{
		return 3;
	}
	cudaStatus = cudaMemcpy(d_line_lat, line_lat, size_double_polyline, cudaMemcpyKind::cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess)
	{
		return 3;
	}
	cudaStatus = cudaMemcpy(d_line_lon, line_lon, size_double_polyline, cudaMemcpyKind::cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess)
	{
		return 3;
	}


	dotarray_near_polyline <<< dimGrid, dimBlock >>> (d_dot_lat, d_dot_lon, d_line_lat, d_line_lon, dot_count, line_count, max_delta, d_dot_result);
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		const char* cudaerr = cudaGetErrorString(cudaStatus);
		//char source[] = (char[])cudaerr;
		int Size;
		while (cudaerr[Size] != '\0') Size++;
		sprintf_s(str, Size, cudaerr);
		return 4;
	}

	cudaStatus = cudaMemcpy(dot_result, d_dot_result, size_int, cudaMemcpyKind::cudaMemcpyDeviceToHost);
	if (cudaStatus != cudaSuccess)
	{
		return 3;
	}

	cudaStatus = cudaFree(d_dot_lat);
	if (cudaStatus != cudaSuccess)
	{
		return 5;
	}
	cudaStatus = cudaFree(d_dot_lon);
	if (cudaStatus != cudaSuccess)
	{
		return 5;
	}	
	cudaStatus = cudaFree(d_line_lat);
	if (cudaStatus != cudaSuccess)
	{
		return 5;
	}	
	cudaStatus = cudaFree(d_line_lon);
	if (cudaStatus != cudaSuccess)
	{
		return 5;
	}	
	cudaStatus = cudaFree(d_dot_result);
	if (cudaStatus != cudaSuccess)
	{
		return 5;
	}

	return 0;
}
