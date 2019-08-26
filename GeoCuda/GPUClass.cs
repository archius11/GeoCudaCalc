using System;
using System.Text;
using System.Runtime.InteropServices;
using System.Collections.Generic;
using System.Linq;
using System.Threading;
using System.Threading.Tasks;

namespace GeoCuda
{
    //------------------------------------ Блок GPU ------------------------------------

    //Инициализация CUDA
    public partial class Calc : IMyClass
    {
        /* return error codes:
        * 0 - success
        * 1 - cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?
        * 2 - cudaMalloc failed!
        * 3 - cudaMemcpy failed!
        * 4 - kernel launch failed!
        * 5 - cudaFree failed!
        * */


        public string errorcode { get; set; }

        [DllImport("cuda_dll.dll", CallingConvention = CallingConvention.Cdecl)]
        private static extern int CudaInitialize();

        void SetErrorCode(int code)
        {
            switch (code)
            {
                case 1: errorcode = "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?"; break;
                case 2: errorcode = "cudaMalloc failed!"; break;
                case 3: errorcode = "cudaMemcpy failed!"; break;
                case 4: errorcode = "kernel launch failed!"; break;
                case 5: errorcode = "cudaFree failed!"; break;
            }
        }

        public bool CudaInit()
        {
            int result = CudaInitialize();
            if (result != 0)
            {
                SetErrorCode(result);
                return false;
            }

            return true;

        }

    }

    //Прямая и обратная геодезические задачи для массива
    public partial class Calc : IMyClass
    {
        [DllImport("cuda_dll.dll", CallingConvention = CallingConvention.Cdecl)]
        private static extern int GetInvertGeo(double[] dot1_lat, double[] dot1_lon, double[] dot2_lat, double[] dot2_lon, double[] dist, double[] azimut, long count);

        [DllImport("cuda_dll.dll", CallingConvention = CallingConvention.Cdecl)]
        private static extern int GetDirectGeo(double[] dot1_lat, double[] dot1_lon, double[] dist, double[] azimut, double[] dot2_lat, double[] dot2_lon, long count);

        public bool GeoInvertProblemArray([MarshalAs(UnmanagedType.SafeArray, SafeArraySubType = VarEnum.VT_R8)] double[] dot1_lat,
                                     [MarshalAs(UnmanagedType.SafeArray, SafeArraySubType = VarEnum.VT_R8)] double[] dot1_lon,
                                     [MarshalAs(UnmanagedType.SafeArray, SafeArraySubType = VarEnum.VT_R8)] double[] dot2_lat,
                                     [MarshalAs(UnmanagedType.SafeArray, SafeArraySubType = VarEnum.VT_R8)] double[] dot2_lon,
                                     long count)
        {
            distance_array = new double[count];
            azimut_array = new double[count];

            int result = GetInvertGeo(dot1_lat, dot1_lon, dot2_lat, dot2_lon, distance_array, azimut_array, count);
            if (result != 0)
            {
                SetErrorCode(result);
                return false;
            }

            for (long i = 0; i < count; i++)
            {
                distance_array[i] = Math.Round(distance_array[i], 2);
                azimut_array[i] = Math.Round(azimut_array[i], 6);
            }

            return true;

        }


        public bool GeoDirectProblemArray([MarshalAs(UnmanagedType.SafeArray, SafeArraySubType = VarEnum.VT_R8)] double[] dot1_lat,
                                     [MarshalAs(UnmanagedType.SafeArray, SafeArraySubType = VarEnum.VT_R8)] double[] dot1_lon,
                                     [MarshalAs(UnmanagedType.SafeArray, SafeArraySubType = VarEnum.VT_R8)] double[] param_dist,
                                     [MarshalAs(UnmanagedType.SafeArray, SafeArraySubType = VarEnum.VT_R8)] double[] param_azimut,
                                     long count)
        {
            dot2_lat_array = new double[count];
            dot2_lon_array = new double[count];


            int result = GetDirectGeo(dot1_lat, dot1_lon, param_dist, param_azimut, dot2_lat_array, dot2_lon_array, count);
            if (result != 0)
            {
                SetErrorCode(result);
                return false;
            }

            for (long i = 0; i < count; i++)
            {
                dot2_lat_array[i] = Math.Round(dot2_lat_array[i], 6);
                dot2_lon_array[i] = Math.Round(dot2_lon_array[i], 6);
            }

            return true;

        }


    }

    //DotArrayCloseToPolylineGPU
    public partial class Calc : IMyClass
    {

        //[DllImport("cuda_dll.dll", CallingConvention = CallingConvention.Cdecl)]
        //private static extern int DotArrayNearPolylineGPU(double[] dot1_lat, double[] dot1_lon, double[] dot2_lat, double[] dot2_lon, double[] dist, double[] azimut, long count);

        //public bool DotArrayCloseToPolylineGPU([MarshalAs(UnmanagedType.SafeArray, SafeArraySubType = VarEnum.VT_R8)] double[] dot_lat,
        //                                    [MarshalAs(UnmanagedType.SafeArray, SafeArraySubType = VarEnum.VT_R8)] double[] dot_lon,
        //                                    [MarshalAs(UnmanagedType.SafeArray, SafeArraySubType = VarEnum.VT_R8)] double[] line_lat,
        //                                    [MarshalAs(UnmanagedType.SafeArray, SafeArraySubType = VarEnum.VT_R8)] double[] line_lon,
        //                                    [MarshalAs(UnmanagedType.SafeArray, SafeArraySubType = VarEnum.VT_INT)] int[] segments,
        //                                    long dot_count,
        //                                    long line_count,
        //                                    float max_delta)
        //{
        //    var sb = new StringBuilder(4096);
        //    IntPtr dot_result = IntPtr.Zero;
        //    int result = DotArrayNearPolylineGPU(dot_lat, dot_lon, line_lat, line_lon, segments, dot_count, line_count, max_delta, ref dot_result, sb);

        //    Marshal.PtrToStructure(dot_result,)
        //    //int result = DotArrayNearPolyline(dot_lat, dot_lon, line_lat, line_lon, dot_count, line_count, max_delta, intarr, sb);
        //    /*if (result != 0)
        //    {
        //        SetErrorCode(result);
        //        errorcode = cudaer;
        //        return false;
        //    }*/
        //    if (result != 0)
        //    {
        //        errorcode = sb.ToString();
        //        return false;
        //    }

        //    for (int i=0; i< dot_count; i++)
        //    {
        //        errorcode = "" + dot_result[0][0] + dot_result[0][1] + dot_result[0][2] + dot_result[1][0] + dot_result[2][1] + dot_result[3][2];
        //    }

        //    return true;
        //}


    }
}
