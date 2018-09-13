using System;
using System.Text;
using System.Runtime.InteropServices;
using System.Collections.Generic;



namespace GeoCuda
{
    [Guid("AC4F4E46-07F8-48C9-BA33-BA8FE0AE75EC")]
    internal interface IMyClass
    {
        [DispId(1)]

        bool GeoInvertProblemArray(double[] dot1_lat, double[] dot1_lon, double[] dot2_lat, double[] dot2_lon, long count);
        bool GeoDirectProblemArray(double[] dot1_lat, double[] dot1_lon, double[] param_dist, double[] param_azimut, long count);
        bool GeoInvertProblem(double dot1_lat, double dot1_lon, double dot2_lat, double dot2_lon);
        bool GeoDirectProblem(double dot1_lat, double dot1_lon, double param_dist, double param_azimut);
        double DistanceToPolyline(double Mlat, double Mlon, double[] line_array_lat, double[] line_array_lon, long n, double dist_limit);
        bool DotCloseToPolyline(double Mlat, double Mlon, double[] line_array_lat, double[] line_array_lon, long n, float max_delta);

        bool DotArrayCloseToPolylineGPU(double[] dot_lat, double[] dot_lon, double[] line_lat, double[] line_lon, long dot_count, long line_count, float max_delta);

        //bool DotArrayCloseToPolylineCPU(double[] dot_lat, double[] dot_lon, double[] line_lat, double[] line_lon, int[] segments, long dot_count, long line_count, float max_delta);

        string[] GetSegmentList();

        bool CudaInit();
        void CleanInvert();
        void CleanDirect();
        double[] distance_array { get; set; }
        double[] azimut_array { get; set; }
        double[] dot2_lat_array { get; set; }
        double[] dot2_lon_array { get; set; }
        double distance { get; set; }
        double azimut { get; set; }
        double lat { get; set; }
        double lon { get; set; }
        string errorcode { get; set; }
        bool[] bool_array { get; set; }

        string[] segment_list { get; set; }


    }

    [Guid("61337B86-DE1F-45D5-A5A8-2B612057E860"), InterfaceType(ComInterfaceType.InterfaceIsIDispatch)]
    public interface IMyEvents
    {
    }

    [Guid("D891EFCF-D6DB-4CB2-943E-A460E6EDDCDB"), ClassInterface(ClassInterfaceType.None), ComSourceInterfaces(typeof(IMyEvents))]
    public class Calc : IMyClass //название нашего класса MyClass
    {
        public double[] distance_array { get; set; }
        public double[] azimut_array { get; set; }
        public double[] dot2_lat_array { get; set; }
        public double[] dot2_lon_array { get; set; }

        public double distance { get; set; }
        public double azimut { get; set; }
        public double lat { get; set; }
        public double lon { get; set; }

        public bool[] bool_array { get; set; }

        public string[] segment_list {get; set; }



        double _PI = 3.14159265358979323846;
        double _PI2 = 1.57079632679489661923;
        double _RAD = 6372795;

        /* return error codes:
         * 0 - success
         * 1 - cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?
         * 2 - cudaMalloc failed!
         * 3 - cudaMemcpy failed!
         * 4 - kernel launch failed!
         * 5 - cudaFree failed!
         * */

        [DllImport("cuda_dll.dll", CallingConvention = CallingConvention.Cdecl)]
        private static extern int GetInvertGeo(double[] dot1_lat, double[] dot1_lon, double[] dot2_lat, double[] dot2_lon, double[] dist, double[] azimut, long count);

        [DllImport("cuda_dll.dll", CallingConvention = CallingConvention.Cdecl)]
        private static extern int GetDirectGeo(double[] dot1_lat, double[] dot1_lon, double[] dist, double[] azimut,  double[] dot2_lat, double[] dot2_lon, long count);

        [DllImport("cuda_dll.dll", CallingConvention = CallingConvention.Cdecl)]
        private static extern int CudaInitialize();

        [DllImport("cuda_dll.dll", CallingConvention = CallingConvention.Cdecl)]
        private static extern int DotArrayNearPolyline(double[] dot_lat, double[] dot_lon, double[] line_lat, double[] line_lon, long dot_count, long line_count, float max_delta, int[] dot_result, StringBuilder str);


        public string errorcode { get; set; }

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

        public void CleanInvert()
        {
            
            distance_array = null;
            azimut_array = null;
            GC.Collect();
            return;
        }
        public void CleanDirect()
        {
            dot2_lat_array = null;
            dot2_lon_array = null;
            GC.Collect();
            return;
        }

        /*[return: MarshalAs(UnmanagedType.SafeArray, SafeArraySubType = VarEnum.VT_BSTR)]
        public  string[] arrtest()
        {
            segment_list = new string[2] { "test", "test2" };

            azimut_array = new double[2] { 1, 2 };



            return segment_list;
        }*/
        [return: MarshalAs(UnmanagedType.SafeArray, SafeArraySubType = VarEnum.VT_BSTR)]
        public string[] GetSegmentList()
        {
            return segment_list;
        }

        public bool GeoInvertProblemArray([MarshalAs(UnmanagedType.SafeArray, SafeArraySubType = VarEnum.VT_R8)] double[] dot1_lat,
                                     [MarshalAs(UnmanagedType.SafeArray, SafeArraySubType = VarEnum.VT_R8)] double[] dot1_lon,
                                     [MarshalAs(UnmanagedType.SafeArray, SafeArraySubType = VarEnum.VT_R8)] double[] dot2_lat,
                                     [MarshalAs(UnmanagedType.SafeArray, SafeArraySubType = VarEnum.VT_R8)] double[] dot2_lon,
                                     long count)
        {
            distance_array = new double[count];
            azimut_array = new double[count];

            int result = GetInvertGeo(dot1_lat, dot1_lon, dot2_lat, dot2_lon, distance_array, azimut_array, count);
            if (result!=0)
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
            if(result!=0)
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

        public bool GeoInvertProblem(double dot1lat,
                                         double dot1lon,
                                         double dot2lat,
                                         double dot2lon)
        {
            _PI = 3.14159265358979323846;
            _PI2 = 1.57079632679489661923;
            _RAD = 6372795;

            double dot1_lat = dot1lat * _PI / 180;  //lat1
            double dot1_lon = dot1lon * _PI / 180;  //lng1
            double dot2_lat = dot2lat * _PI / 180;  //lat2
            double dot2_lon = dot2lon * _PI / 180;  //lng2

            double cl1, cl2, sl1, sl2, delta, cdelta, sdelta;

            cl1 = Math.Cos(dot1_lat);
            cl2 = Math.Cos(dot2_lat);
            sl1 = Math.Sin(dot1_lat);
            sl2 = Math.Sin(dot2_lat);
            delta = dot2_lon - dot1_lon;
            cdelta = Math.Cos(delta);
            sdelta = Math.Sin(delta);

            double x, y, z, ad, z2;
            y = Math.Sqrt(Math.Pow(cl2 * sdelta, 2) + Math.Pow(cl1 * sl2 - sl1 * cl2 * cdelta, 2));
            x = sl1 * sl2 + cl1 * cl2 * cdelta;
            ad = Math.Atan(y/x);

            distance = ad * _RAD; //current distance

            x = (cl1 * sl2) - (sl1 * cl2 * cdelta);
            y = sdelta * cl2;

            z = 0; //
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
                z = Math.Atan(-y/x) * 180 / _PI;
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
            anglerad2 = z2 - ((2 * _PI) * Math.Floor(z2 / (2 * _PI)));
            azimut = anglerad2 * 180 / _PI;

            distance = Math.Round(distance, 2);
            azimut = Math.Round(azimut, 6);
            return true;
        }

        public bool GeoDirectProblem(double dot1lat,
                                         double dot1lon,
                                         double paramdist,
                                         double paramazimut)
        {

            _PI = 3.14159265358979323846;
            _PI2 = 1.57079632679489661923;
            _RAD = 6372795;

            double dot1_lat = dot1lat * _PI / 180;
            double dot1_lon = dot1lon * _PI / 180;
            double param_dist = paramdist / _RAD;
            double param_azimut = paramazimut * _PI / 180;

            double[] pt2 = new double[2];
            double[] pt = new double[] { dot1_lat, dot1_lon};
            SphereDirect(pt, param_azimut, param_dist, pt2);

            pt2[0] = pt2[0] * 180 / _PI;
            pt2[1] = pt2[1] * 180 / _PI;

            if (pt2[0] < 0)
                pt2[0] += 180;

            if (pt2[1] < 0)
                pt2[1] += 180;

            lat = pt2[0];
            lon = pt2[1];

            lat = Math.Round(lat, 6);
            lon = Math.Round(lon, 6);

            return true;

        }

        public bool DotArrayCloseToPolylineGPU([MarshalAs(UnmanagedType.SafeArray, SafeArraySubType = VarEnum.VT_R8)] double[] dot_lat,
                                            [MarshalAs(UnmanagedType.SafeArray, SafeArraySubType = VarEnum.VT_R8)] double[] dot_lon,
                                            [MarshalAs(UnmanagedType.SafeArray, SafeArraySubType = VarEnum.VT_R8)] double[] line_lat,
                                            [MarshalAs(UnmanagedType.SafeArray, SafeArraySubType = VarEnum.VT_R8)] double[] line_lon,
                                            long dot_count,
                                            long line_count,
                                            float max_delta)
        {
            bool_array = new bool[dot_count];
            int[] intarr = new int[dot_count];

            var sb = new StringBuilder(4096);
            int result = DotArrayNearPolyline(dot_lat, dot_lon, line_lat, line_lon, dot_count, line_count, max_delta, intarr, sb);
            /*if (result != 0)
            {
                SetErrorCode(result);
                errorcode = cudaer;
                return false;
            }*/
            if (result != 0)
            {
                errorcode = sb.ToString();
                return false;
            }

            for (int i=0; i< dot_count; i++)
            {
                bool_array[i] = intarr[i] == 1;
            }

            return true;
        }

        public bool CudaInit()
        {
            int result = CudaInitialize();
            if (result!=0)
            {
                SetErrorCode(result);
                return false;
            }

            return true;

        }

        double hypot(double x, double y)
        {
            return Math.Sqrt(x * x + y * y);
        }

        void CartToSpher(double[] x, double[] y)
        {
            double p;

            p = hypot(x[0], x[1]);
            y[1] = Math.Atan(x[1]/x[0]);
            y[0] = Math.Atan(x[2]/p);

            return;
        }


        void SpherToCart(double[] y, double[] x)
        {
            double p;

            p = Math.Cos(y[0]);
            x[2] = Math.Sin(y[0]);
            x[1] = p * Math.Sin(y[1]);
            x[0] = p * Math.Cos(y[1]);

            return;
        }


        void Rotate2(double[] x, double a)
        {
            double c, s, xj;

            c = Math.Cos(a);
            s = Math.Sin(a);
            xj = x[0] * c + x[1] * s;
            x[1] = -x[0] * s + x[1] * c;
            x[0] = xj;

            return;
        }

        void Rotate1(double[] x, double a)
        {
            double c, s, xj;

            c = Math.Cos(a);
            s = Math.Sin(a);
            xj = x[2] * c + x[0] * s;
            x[0] = -x[2] * s + x[0] * c;
            x[2] = xj;

            return;
        }


        void SphereDirect(double[] pt1, double azi, double dist, double[] pt2)
        {
            double[] pt = new double[2];
            double[] x = new double[3];

            pt[0] = _PI2 - dist;
            pt[1] = _PI - azi;

            SpherToCart(pt, x);               // сферические -> декартовы
            Rotate1(x, pt1[0] - _PI2); // первое вращение
            Rotate2(x, -pt1[1]);           // второе вращение
            CartToSpher(x, pt2);           // декартовы -> сферические 

        }

        void SpherToCartR(double[] y, double[] x)
        {
            double y_lat = y[0] * _PI / 180;
            double y_lon = y[1] * _PI / 180;

            double cos_y_lat = Math.Cos(y_lat);

            x[0] = _RAD * cos_y_lat * Math.Cos(y_lon);
            x[1] = _RAD * cos_y_lat * Math.Sin(y_lon);
            x[2] = _RAD * Math.Sin(y_lat);

        }

        public bool DotCloseToPolyline(double Mlat, double Mlon,
                                        [MarshalAs(UnmanagedType.SafeArray, SafeArraySubType = VarEnum.VT_R8)] double[] line_array_lat,
                                        [MarshalAs(UnmanagedType.SafeArray, SafeArraySubType = VarEnum.VT_R8)] double[] line_array_lon,
                                        long n, float max_delta)
        {
            _PI = 3.14159265358979323846;
            _PI2 = 1.57079632679489661923;
            _RAD = 6372795;

            //long n = line_array_lat.Length;

            //dist_limit = -1;
            double[] m_sph = new double[2] { Mlat, Mlon };
            double[] m_dec = new double[3];

            SpherToCartR(m_sph, m_dec);

            double disttoline;
            double d_max_delta = max_delta;
            bool dotisclose = false;

            for (long i = 0; i <= n - 2; i++)
            {
                //line dots
                double[] a_sph = new double[2] { line_array_lat[i], line_array_lon[i] };
                double[] b_sph = new double[2] { line_array_lat[i + 1], line_array_lon[i + 1] };

                if (LineIsTooFar(m_sph, a_sph, b_sph, d_max_delta))
                {
                    continue;
                }

                double[] a_dec = new double[3];
                double[] b_dec = new double[3];

                SpherToCartR(a_sph, a_dec);
                SpherToCartR(b_sph, b_dec);

                disttoline = DistanceToLine(m_dec, a_dec, b_dec, i == n - 2);

                if (disttoline < max_delta)
                {
                    dotisclose = true;
                    break;
                }
            }

            return dotisclose;
        }



        /*
        {"#",51e7a0d2-530b-11d4-b98a-008048da3034,
        {4,
        {"N",1},
        {"N",4024},
        {"N",1},
        {"N",4024}
        }
        }

        {"#",51e7a0d2-530b-11d4-b98a-008048da3034,
        {0}
        }

        */

        int[] DotCloseToSegments(double Mlat, double Mlon,
                                double[][][] segments_array, float max_delta)
        {
            _PI = 3.14159265358979323846;
            _PI2 = 1.57079632679489661923;
            _RAD = 6372795;

            //long n = line_array_lat.Length;

            //dist_limit = -1;
            double[] m_sph = new double[2] { Mlat, Mlon };
            double[] m_dec = new double[3];

            SpherToCartR(m_sph, m_dec);

            double disttoline;
            double d_max_delta = max_delta;
            List<int> segment_list = new List<int>() { };
            int segments_count = segments_array.Length;

            for (int l = 0; l<segments_count; l++)
            {
                double[] segment_lats = segments_array[l][0];
                double[] segment_lons = segments_array[l][1];

                int n = segment_lats.Length;
                //bool segmentclose = false;
                for (long i = 0; i <= n - 2; i++)
                {
                    //line dots
                    double[] a_sph = new double[2] { segment_lats[i], segment_lons[i] };
                    double[] b_sph = new double[2] { segment_lats[i + 1], segment_lons[i + 1] };

                    if (LineIsTooFar(m_sph, a_sph, b_sph, d_max_delta))
                    {
                        continue;
                    }

                    double[] a_dec = new double[3];
                    double[] b_dec = new double[3];

                    SpherToCartR(a_sph, a_dec);
                    SpherToCartR(b_sph, b_dec);

                    disttoline = DistanceToLine(m_dec, a_dec, b_dec, i == n - 2);

                    if (disttoline < max_delta)
                    {
                        segment_list.Add(l);
                        break;
                    }
                }
            }

            return segment_list.ToArray();
        }

        public bool DotArrayCloseToPolylineCPU([MarshalAs(UnmanagedType.SafeArray, SafeArraySubType = VarEnum.VT_R8)] double[] dot_lat,
                                            [MarshalAs(UnmanagedType.SafeArray, SafeArraySubType = VarEnum.VT_R8)] double[] dot_lon,
                                            [MarshalAs(UnmanagedType.SafeArray, SafeArraySubType = VarEnum.VT_R8)] double[] line_lat,
                                            [MarshalAs(UnmanagedType.SafeArray, SafeArraySubType = VarEnum.VT_R8)] double[] line_lon,
                                            [MarshalAs(UnmanagedType.SafeArray, SafeArraySubType = VarEnum.VT_INT)] int[] segments,
                                            long dot_count,
                                            long line_count,
                                            float max_delta)
        {
            List<List<double>[]> segments_list = new List<List<double>[] > { }; //о_О список массивов списков
            //segments_arrays.Add(new List<int>[2]); //добавили новый сегмент (массив из 2 списков ширин и долгот)
            //segments_arrays[0][0].Add(1);     //первый элемент списка, первый список массива, добавили в список
            

            //int[][][] segments_arrays; //массив сегментов, каждый сегмент это массив из 2 элементов - массивов, массив ширина и массив долгота
            //segments_arrays  = new int[1][][]; //первый сегмент
            //segments_arrays[0] = new int[2][]; //инициализация массива из 2 массивов широт долгот
            //segments_arrays[0][0] = new int[1]; //массив широт, первая широта
            //segments_arrays[0][1] = new int[1]; //массив широт, первая широта

            int tek_seg = -1; //текущее имя сегмента, -1 неинициализирован
            int segments_count = 0; //текущее количество сегментов
            int dots_in_segment_count = 0; //количество точек в текущем сегменте

            for (int i = 0; i < line_count; i++) //цикл по каждой точке трека (количество элементов линии и сегментов одинаково)
            {
                //если сменился номер сегмента, то добавим новый сегмент
                if (segments[i]!=tek_seg)
                {
                    tek_seg = segments[i];
                    segments_list.Add(new List<double>[2]); //добавили новый сегмент в виде массива из 2 списков чисел
                    segments_list[segments_count][0] = new List<double>() { };
                    segments_list[segments_count][1] = new List<double>() { };
                    segments_count += 1; //+1 сегмент
                    dots_in_segment_count = 0; //новый сегмент без точек
                }
                segments_list[segments_count - 1][0].Add(line_lat[i]); //в текущий список добавили широту
                segments_list[segments_count - 1][1].Add(line_lon[i]); //в текущий список добавили долготу

                dots_in_segment_count += 1;

            }


            double[][][] segments_array;
            segments_array = new double[segments_count][][];
            for (int i = 0; i<segments_count; i++)
            {
                segments_array[i] = new double[2][];
                segments_array[i][0] = segments_list[i][0].ToArray(); //массив широт
                segments_array[i][1] = segments_list[i][1].ToArray(); //массив долгот
            }



            segment_list = new string[dot_count];

            int[] closesegments;
            string segments_array_string;
            for (int i=0; i<dot_count; i++)
            {
                closesegments = DotCloseToSegments(dot_lat[i], dot_lon[i], segments_array, max_delta); //int[]
                if(closesegments.Length==0)
                {
                    segments_array_string = "{\"#\",51e7a0d2-530b-11d4-b98a-008048da3034,{0}}";
                }
                else
                {
                    /*{
                        "#",51e7a0d2- 530b - 11d4 - b98a - 008048da3034,
                        {
                                            4,
                        { "N",1},
                        { "N",4024},
                        { "N",1},
                        { "N",4024}
                        }
                    }*/
                    int closesegments_Length = closesegments.Length;
                    segments_array_string = "{\"#\",51e7a0d2-530b-11d4-b98a-008048da3034,{" + closesegments_Length + ",";
                    for(int segmentnum = 0; segmentnum < closesegments_Length; segmentnum++)
                    {
                        segments_array_string += "{\"N\"," + closesegments[segmentnum] + "}";
                        if(segmentnum<closesegments_Length-1)
                        {
                            segments_array_string += ",";

                        }
                    }
                    segments_array_string += "}}";
                }

                segment_list[i] = segments_array_string;

            }


            return true;
        }

        public double DistanceToPolyline(double Mlat, double Mlon,
                                        [MarshalAs(UnmanagedType.SafeArray, SafeArraySubType = VarEnum.VT_R8)] double[] line_array_lat,
                                        [MarshalAs(UnmanagedType.SafeArray, SafeArraySubType = VarEnum.VT_R8)] double[] line_array_lon,
                                        long n, double dist_limit)
        {
            _PI = 3.14159265358979323846;
            _PI2 = 1.57079632679489661923;
            _RAD = 6372795;

            //long n = line_array_lat.Length;

            //dist_limit = -1;
            double[] m_sph = new double[2] { Mlat, Mlon };
            double[] m_dec = new double[3];

            SpherToCartR(m_sph, m_dec);

            double disttoline, disttopolyline;
            disttopolyline = -1;

            for (long i=0; i<= n - 2; i++)
            {
                //line dots
                double[] a_sph = new double[2] { line_array_lat[i], line_array_lon[i] };
                double[] b_sph = new double[2] { line_array_lat[i+1], line_array_lon[i+1] };

                if (dist_limit != -1)
                {
                    if(LineIsTooFar(m_sph, a_sph, b_sph, dist_limit))
                    {
                        continue;
                    }
                }

                double[] a_dec = new double[3];
                double[] b_dec = new double[3];

                SpherToCartR(a_sph, a_dec);
                SpherToCartR(b_sph, b_dec);

                disttoline = DistanceToLine(m_dec, a_dec, b_dec, i== n - 2);

                if(disttopolyline==-1)
                {
                    disttopolyline = disttoline;
                }
                else
                {
                    disttopolyline = Math.Min(disttopolyline, disttoline);
                }
            }

            return disttopolyline;
        }

        bool LineIsTooFar(double[] M, double[] A, double[] B, double limit)
        {
            double lat_lag = 0.00001 * limit * 1.2;
            double lon_lag = 0.00003 * limit * 1.2;

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

        double DistanceToLine(double[] m_dec, double[] a_dec, double[] b_dec, bool calc_mb)
        {
            if (a_dec[0] == b_dec[0] &&
                a_dec[1] == b_dec[1] &&
                a_dec[2] == b_dec[2])
                {
                    LineLength(m_dec, a_dec);
                }

            double plane_A, plane_B, plane_C;
            plane_A = a_dec[1] * b_dec[2] - a_dec[2] * b_dec[1];
            plane_B = a_dec[2] * b_dec[0] - a_dec[0] * b_dec[2];
            plane_C = a_dec[0] * b_dec[1] - a_dec[1] * b_dec[0];

            double d, MK_length, MA_length, MB_length, minlength;

            if (M_ProjectionOnPlane(m_dec, plane_A, plane_B, plane_C, a_dec, b_dec))
            {
                d = Math.Abs(plane_A*m_dec[0] + plane_B * m_dec[1] + plane_C*m_dec[2]) / Math.Sqrt(Math.Pow(plane_A,2) + Math.Pow(plane_B, 2) + Math.Pow(plane_C, 2));
                MK_length = _RAD * Math.Asin(d / _RAD);
                minlength = MK_length;
            }
            else
            {
                MA_length = LineLength(m_dec, a_dec);
                if(calc_mb)
                {
                    MB_length = LineLength(m_dec, b_dec);
                    minlength = Math.Min(MA_length, MB_length);
                }
                else
                {
                    minlength = MA_length;
                }
            }

            return minlength;
        }

        double LineLength(double[] a, double[] b)
        {
            double K1, K2, K3, K4;

            K1 = b[0] * a[0] + b[1] * a[1] + b[2] * a[2];
            K2 = Math.Sqrt(Math.Pow(b[0], 2) + Math.Pow(b[1], 2) + Math.Pow(b[2], 2));
            K3 = Math.Sqrt(Math.Pow(a[0], 2) + Math.Pow(a[1], 2) + Math.Pow(a[2], 2));
            K4 = K1 / (K2 * K3);

            if (K4 > 1) K4 = 1;

            return _RAD * Math.Acos(K4);
        }

        bool M_ProjectionOnPlane(double[] m_dec, double plane_A, double plane_B, double plane_C, double[] a_dec, double[] b_dec)
        {
            double t = -(plane_A * m_dec[0] + plane_B * m_dec[1] + plane_C * m_dec[2]) / (Math.Pow(plane_A, 2) + Math.Pow(plane_B, 2) + Math.Pow(plane_C, 2));

            double[] dot_on_plane = new double[3];
            dot_on_plane[0] = plane_A * t + m_dec[0];
            dot_on_plane[1] = plane_B * t + m_dec[1];
            dot_on_plane[2] = plane_C * t + m_dec[2];

            double[] dot_k = new double[3];
            double K = Math.Sqrt(Math.Pow(dot_on_plane[0], 2) + Math.Pow(dot_on_plane[1], 2) + Math.Pow(dot_on_plane[2], 2));
            dot_k[0] = (_RAD * dot_on_plane[0]) / K;
            dot_k[1] = (_RAD * dot_on_plane[1]) / K;
            dot_k[2] = (_RAD * dot_on_plane[2]) / K;

            double line_AB = LineLength(a_dec, b_dec);
            double line_AK = LineLength(a_dec, dot_k);
            double line_BK = LineLength(b_dec, dot_k);

            return Math.Abs(line_AK + line_BK - line_AB) < 0.01;
        }


    }
}
