using System;
using System.Text;
using System.Runtime.InteropServices;
using System.Collections.Generic;
using System.Linq;
using System.Threading;
using System.Threading.Tasks;
using System.Globalization;


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
        bool DotInsidePolygons(double DotLat, double DotLon, double[] PolyLat, double[] PolyLon);
        bool DotsInsidePolygon(double[] DotLats, double[] DotLons, double[] PolyLat, double[] PolyLon);
        bool Test2dArrayMashalling(double[,] arr);
        bool TestArrayOfArraysMashalling(Array arr);
        string GetDoubleArrayAsString(double[] doubles);

        //bool DotArrayCloseToPolylineGPU(double[] dot_lat, double[] dot_lon, double[] line_lat, double[] line_lon, int[] segments, long dot_count, long line_count, float max_delta);

        bool DotArrayCloseToPolylineCPU(double[] dot_lat, double[] dot_lon, double[] line_lat, double[] line_lon, int[] segments, long dot_count, long line_count, float max_delta);

        string[] GetSegmentList();

        double[,] TestGet2dArray();

        string[] TestGetArrayOfArrays(Array arr);
        


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
        bool IsInside { get; set; }
        bool[] IsInsideArray { get; set; }

        int[,] TrackIsInsideArray { get; set; }

        string[] segment_list { get; set; }

        double[][] test_return_arrayofarrays { get; set; }

    }

    [Guid("61337B86-DE1F-45D5-A5A8-2B612057E860"), InterfaceType(ComInterfaceType.InterfaceIsIDispatch)]
    public interface IMyEvents
    {
    }

    //массивы переменные класса
    [Guid("D891EFCF-D6DB-4CB2-943E-A460E6EDDCDB"), ClassInterface(ClassInterfaceType.None), ComSourceInterfaces(typeof(IMyEvents))]
    public partial class Calc : IMyClass
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

        public bool IsInside { get; set; }
        public bool[] IsInsideArray { get; set; }
        public int[,] TrackIsInsideArray { get; set; }

        public double[][] test_return_arrayofarrays { get; set; }


        double _PI = 3.14159265358979323846;
        double _PI2 = 1.57079632679489661923;
        double _RAD = 6372795;



    }


    //------------------------------------ Блок CPU ------------------------------------

    //Универсальные функции
    public partial class Calc : IMyClass
    {


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

        double hypot(double x, double y)
        {
            return Math.Sqrt(x * x + y * y);
        }
        
        void CartToSpher(double[] x, double[] y)
        {
            double p;

            p = hypot(x[0], x[1]);
            y[1] = Math.Atan(x[1] / x[0]);
            y[0] = Math.Atan(x[2] / p);

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

        double DistanceToLine(double[] m_dec, double[] a_dec, double[] b_dec, bool calc_mb)
        {
            if (a_dec[0] == b_dec[0] && //TODO: V3024 https://www.viva64.com/en/w/v3024/ An odd precise comparison: a_dec[0] == b_dec[0]. Consider using a comparison with defined precision: Math.Abs(A - B) < Epsilon.
                a_dec[1] == b_dec[1] && //TODO: V3024 https://www.viva64.com/en/w/v3024/ An odd precise comparison: a_dec[1] == b_dec[1]. Consider using a comparison with defined precision: Math.Abs(A - B) < Epsilon.
                a_dec[2] == b_dec[2]) //TODO: V3024 https://www.viva64.com/en/w/v3024/ An odd precise comparison: a_dec[2] == b_dec[2]. Consider using a comparison with defined precision: Math.Abs(A - B) < Epsilon.
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
                d = Math.Abs(plane_A * m_dec[0] + plane_B * m_dec[1] + plane_C * m_dec[2]) / Math.Sqrt(Math.Pow(plane_A, 2) + Math.Pow(plane_B, 2) + Math.Pow(plane_C, 2));
                MK_length = _RAD * Math.Asin(d / _RAD);
                minlength = MK_length;
            }
            else
            {
                MA_length = LineLength(m_dec, a_dec);
                if (calc_mb)
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

        bool DoublesAreEqual(double A, double B)
        {
            //return Math.Abs(A - B) < 0.000001;
            return A == B;

        }

        bool DoublesAreLessOrEqual(double A, double B)
        {
            //return Math.Abs(A - B) < 0.000001 || A < B;
            return A <= B;

        }

        bool DoublesAreMoreOrEqual(double A, double B)
        {
            //return Math.Abs(A - B) < 0.000001 || A > B;
            return A >= B;

        }

        public string GetDoubleArrayAsString(double[] doubles)
        {
            string ret_string;
            string specifier;
            specifier = "G";

            if (doubles.Length == 0)
            {
                ret_string = "{\"#\",51e7a0d2-530b-11d4-b98a-008048da3034,{0}}";
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
                int arr_length = doubles.Length;
                ret_string = "{\"#\",51e7a0d2-530b-11d4-b98a-008048da3034,{" + arr_length + ",";
                for (int i = 0; i < arr_length; i++)
                {
                    ret_string += "{\"N\"," + doubles[i].ToString(specifier, CultureInfo.InvariantCulture) + "}";
                    if (i < arr_length - 1)
                    {
                        ret_string += ",";

                    }
                }
                ret_string += "}}";
            }

            return ret_string;

        }

        public string GetIntArrayAsString(int[] ints)
        {
            string ret_string;

            if (ints.Length == 0)
            {
                ret_string = "{\"#\",51e7a0d2-530b-11d4-b98a-008048da3034,{0}}";
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
                int arr_length = ints.Length;
                ret_string = "{\"#\",51e7a0d2-530b-11d4-b98a-008048da3034,{" + arr_length + ",";
                for (int i = 0; i < arr_length; i++)
                {
                    ret_string += "{\"N\"," + ints[i] + "}";
                    if (i < arr_length - 1)
                    {
                        ret_string += ",";

                    }
                }
                ret_string += "}}";
            }

            return ret_string;

        }
    }

    //прямая и обратная геодезические задачи
    public partial class Calc : IMyClass
    {
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
            ad = Math.Atan(y / x);

            distance = ad * _RAD; //current distance

            x = (cl1 * sl2) - (sl1 * cl2 * cdelta);
            y = sdelta * cl2;

            z = 0; //
            if (x == 0) //TODO: V3024 https://www.viva64.com/en/w/v3024/ An odd precise comparison: x == 0. Consider using a comparison with defined precision: Math.Abs(A - B) < Epsilon.
            {
                if (y > 0)
                    z = -90;
                else if (y < 0)
                    z = 90;
                else if (y == 0) //TODO: V3024 https://www.viva64.com/en/w/v3024/ An odd precise comparison: y == 0. Consider using a comparison with defined precision: Math.Abs(A - B) < Epsilon.
                    z = 0;
            }
            else
            {
                z = Math.Atan(-y / x) * 180 / _PI;
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
            double[] pt = new double[] { dot1_lat, dot1_lon };
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

    }

    //Точка близка к полилинии
    public partial class Calc : IMyClass
    {
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

    }

    //DotArrayCloseToPolylineCPU
    public partial class Calc : IMyClass
    {
        //массив 1с
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

        [return: MarshalAs(UnmanagedType.SafeArray, SafeArraySubType = VarEnum.VT_BSTR)]
        public string[] GetSegmentList()
        {
            return segment_list;
        }

        public bool DotArrayCloseToPolylineCPU([MarshalAs(UnmanagedType.SafeArray, SafeArraySubType = VarEnum.VT_R8)] double[] dot_lat, //массив широта
                                            [MarshalAs(UnmanagedType.SafeArray, SafeArraySubType = VarEnum.VT_R8)] double[] dot_lon,    //массив долгота
                                            [MarshalAs(UnmanagedType.SafeArray, SafeArraySubType = VarEnum.VT_R8)] double[] line_lat,   //широта полилинии
                                            [MarshalAs(UnmanagedType.SafeArray, SafeArraySubType = VarEnum.VT_R8)] double[] line_lon,   //долгота полилинии
                                            [MarshalAs(UnmanagedType.SafeArray, SafeArraySubType = VarEnum.VT_INT)] int[] segments,     //номера сегментов
                                            long dot_count,                                                                             //количество точек
                                            long line_count,                                                                            //количество точек в полилинии
                                            float max_delta)                                                                            //макс погрешность близоты
        {

            List<List<double>[]> segments_list = new List<List<double>[]> { }; //о_О список массивов списков
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
            List<int> segments_names = new List<int>() { };



            for (int i = 0; i < line_count; i++) //цикл по каждой точке трека (количество элементов линии и сегментов одинаково)
            {
                //если сменился номер сегмента, то добавим новый сегмент
                if (segments[i] != tek_seg)
                {
                    segments_names.Add(segments[i]);
                    tek_seg = segments[i];
                    segments_list.Add(new List<double>[2]); //добавили новый сегмент в виде массива из 2 списков чисел и список из 4 максмины массива
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
            for (int i = 0; i < segments_count; i++)
            {
                segments_array[i] = new double[3][];
                segments_array[i][0] = segments_list[i][0].ToArray(); //массив широт
                segments_array[i][1] = segments_list[i][1].ToArray(); //массив долгот
                segments_array[i][2] = new double[4];


                segments_array[i][2][0] = segments_array[i][0].Max() + (0.00001 * max_delta * 1.2);
                segments_array[i][2][1] = segments_array[i][0].Min() - (0.00001 * max_delta * 1.2);
                segments_array[i][2][2] = segments_array[i][1].Max() + (0.00003 * max_delta * 1.2);
                segments_array[i][2][3] = segments_array[i][1].Min() - (0.00003 * max_delta * 1.2);
            }


            segment_list = new string[dot_count];

            int[] closesegments;
            string segments_array_string;
            for (int i = 0; i < dot_count; i++)
            {
                closesegments = DotCloseToSegments(dot_lat[i], dot_lon[i], segments_array, max_delta, segments_names); //int[]
                if (closesegments.Length == 0)
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
                    for (int segmentnum = 0; segmentnum < closesegments_Length; segmentnum++)
                    {
                        segments_array_string += "{\"N\"," + closesegments[segmentnum] + "}";
                        if (segmentnum < closesegments_Length - 1)
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

            for (long i = 0; i <= n - 2; i++)
            {
                //line dots
                double[] a_sph = new double[2] { line_array_lat[i], line_array_lon[i] };
                double[] b_sph = new double[2] { line_array_lat[i + 1], line_array_lon[i + 1] };

                if (dist_limit != -1) //TODO: V3024 https://www.viva64.com/en/w/v3024/ An odd precise comparison: dist_limit != -1. Consider using a comparison with defined precision: Math.Abs(A - B) > Epsilon.
                {
                    if (LineIsTooFar(m_sph, a_sph, b_sph, dist_limit))
                    {
                        continue;
                    }
                }

                double[] a_dec = new double[3];
                double[] b_dec = new double[3];

                SpherToCartR(a_sph, a_dec);
                SpherToCartR(b_sph, b_dec);

                disttoline = DistanceToLine(m_dec, a_dec, b_dec, i == n - 2);

                if (disttopolyline == -1) //TODO: V3024 https://www.viva64.com/en/w/v3024/ An odd precise comparison: disttopolyline==-1. Consider using a comparison with defined precision: Math.Abs(A - B) < Epsilon.
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

        int[] DotCloseToSegments(double Mlat, double Mlon,
                                double[][][] segments_array, float max_delta, List<int> segments_names)
        {
            //_PI = 3.14159265358979323846;
            //_PI2 = 1.57079632679489661923;
            //_RAD = 6372795;

            //long n = line_array_lat.Length;

            //dist_limit = -1;
            double[] m_sph = new double[2] { Mlat, Mlon };
            double[] m_dec = new double[3];

            SpherToCartR(m_sph, m_dec);

            double disttoline;
            double d_max_delta = max_delta;
            List<int> segment_list = new List<int>() { };
            int segments_count = segments_array.Length;


            ParallelOptions parOps = new ParallelOptions();
            parOps.MaxDegreeOfParallelism = Environment.ProcessorCount;


            Parallel.For(0, segments_count, parOps, l =>
            {
                int segment_name = segments_names[l];
                double[] segment_lats = segments_array[l][0];
                double[] segment_lons = segments_array[l][1];

                double maxlat, minlat, maxlon, minlon;

                maxlat = segments_array[l][2][0];
                minlat = segments_array[l][2][1];
                maxlon = segments_array[l][2][2];
                minlon = segments_array[l][2][3];

                if (!(Mlat >= minlat &&
                    Mlat <= maxlat &&
                    Mlon >= minlon &&
                    Mlon <= maxlon))
                {
                    return;
                }

                int n = segment_lats.Length;

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
                        segment_list.Add(segment_name);
                        break;
                    }
                }


            });

            return segment_list.ToArray();
        }
    }

    //точка находится внутри полилинии
    //точка внутри геозоны. Геозона - массив широт и долгот
    public partial class Calc : IMyClass
    {

        public bool DotInsidePolygons(double DotLat, double DotLon,
                                            [MarshalAs(UnmanagedType.SafeArray, SafeArraySubType = VarEnum.VT_R8)] double[] PolyLat, //массив широта
                                            [MarshalAs(UnmanagedType.SafeArray, SafeArraySubType = VarEnum.VT_R8)] double[] PolyLon)  //массив долгота)
        {
            IsInside = DotInsidePolygon(DotLat, DotLon, PolyLat, PolyLon);

            return true;
        }

        public bool DotsInsidePolygon([MarshalAs(UnmanagedType.SafeArray, SafeArraySubType = VarEnum.VT_R8)] double[] DotLats,
                                       [MarshalAs(UnmanagedType.SafeArray, SafeArraySubType = VarEnum.VT_R8)] double[] DotLons,
                                            [MarshalAs(UnmanagedType.SafeArray, SafeArraySubType = VarEnum.VT_R8)] double[] PolyLat, //массив широта
                                            [MarshalAs(UnmanagedType.SafeArray, SafeArraySubType = VarEnum.VT_R8)] double[] PolyLon)  //массив долгота)
        {
            int count = DotLats.Length;

            if (count == 0) return false;

            List<int[]> Phases = new List<int[]>();

            bool DotInside = false;
            bool newDotInside;

            for (int i=0; i < count; i++)
            {
                newDotInside = DotInsidePolygon(DotLats[i], DotLons[i], PolyLat, PolyLon);

                if (DotInside != newDotInside)
                {
                    int[] newRow = new int[2];
                    newRow[0] = i; //phase index
                    newRow[1] = newDotInside ? 1 : 0; //1:in, 0:out
                    Phases.Add(newRow);
                }

                DotInside = newDotInside;
            }

            TrackIsInsideArray = new int[Phases.Count, 2];
            for (int i=0; i<Phases.Count; i++)
            {
                TrackIsInsideArray[i, 0] = Phases[i][0];
                TrackIsInsideArray[i, 1] = Phases[i][1];
            }

            return true;
        }


        bool DotInsidePolygon(double mLat, double mLon, double[] lat_array, double[] lon_array)
        {
            double LatMin = lat_array.Min();
            double LatMax = lat_array.Max();
            double LonMin = lon_array.Min();
            double LonMax = lon_array.Max();

            if (!(DoublesAreMoreOrEqual(mLat, LatMin) &&
                DoublesAreLessOrEqual(mLat, LatMax) &&
                DoublesAreMoreOrEqual(mLon, LonMin) &&
                DoublesAreLessOrEqual(mLon, LonMax)))
                return false;

            int IntersectionRight = 0;
            int IntersectionLeft = 0;

            int array_dot_count = lat_array.Length;

            double[] begin_dot = new double[2];
            double[] last_dot = new double[2];

            for (int i = 0; i < array_dot_count; i++)
            {
                begin_dot[0] = lat_array[i]; //0 - широта
                begin_dot[1] = lon_array[i]; //1 - долгота

                if (i != 0)
                {
                    if (((DoublesAreLessOrEqual(begin_dot[1], mLon)) && (last_dot[1] > mLon)) ||
                            ((begin_dot[1] > mLon) && (DoublesAreLessOrEqual(last_dot[1], mLon))))
                    {
                        double IntersectionLat = begin_dot[0] + ((last_dot[0] - begin_dot[0]) / (last_dot[1] - begin_dot[1])) * (mLon - begin_dot[1]);
                        if (IntersectionLat < mLat)
                        {
                            IntersectionLeft++;
                        }
                        if (IntersectionLat > mLat)
                        {
                            IntersectionRight++;
                        }
                    }
                }
                last_dot[0] = begin_dot[0];
                last_dot[1] = begin_dot[1];
            }

            begin_dot = new double[] { lat_array[0], lon_array[0] };
            last_dot = new double[] { lat_array[array_dot_count - 1], lon_array[array_dot_count - 1] };

            if ((begin_dot[0] != last_dot[0]) || (begin_dot[1] != last_dot[1]))
            {
                if (((DoublesAreLessOrEqual(last_dot[1], mLon)) && (begin_dot[1] > mLon)) ||
                    ((last_dot[1] > mLon) && (DoublesAreLessOrEqual(begin_dot[1], mLon))))
                {
                    double IntersectionLat = last_dot[0] + ((begin_dot[0] - last_dot[0]) / (begin_dot[1] - last_dot[1])) * (mLon - last_dot[1]);
                    if (IntersectionLat < mLat)
                        IntersectionLeft++;
                    if (IntersectionLat > mLat)
                        IntersectionRight++;
                }
            }

            return ((IntersectionLeft % 2 == 1) && (IntersectionRight % 2 == 1));
        }
    }

    public partial class Calc : IMyClass
    {

        public bool Test2dArrayMashalling([MarshalAs(UnmanagedType.SafeArray, SafeArraySubType = VarEnum.VT_R8)] double[,] arr)
        {
            int cols = arr.GetLength(0); //columns
            int rows = arr.GetLength(1); //rows
            errorcode = "" + rows + ":" + cols + ": ";
            for (int row = 0; row < rows; row++)
                for (int col = 0; col < cols; col++)
                    errorcode += arr[col, row] + ", ";

            //2:3: 1, 2, 3, 4, 5, 6,   

            return true;
        }

        public bool TestArrayOfArraysMashalling([MarshalAs(UnmanagedType.SafeArray, SafeArraySubType = VarEnum.VT_VARIANT)] Array arr)
        {
            int arr_num = arr.Length; //2 вложенных массива
            errorcode = "good ";

            //double[][] main_Array = new double[arr_num][]; //массив из 2 массивов
            test_return_arrayofarrays = new double[arr_num][];

            //int cnt = 0;
            for (int i=0; i< arr_num; i++)  //по каждому массиву
            {
                double[] d1 = (double[]) arr.GetValue(i); //получение i-массива
                test_return_arrayofarrays[i] = new double[d1.Length]; 
                Array.Copy(d1, test_return_arrayofarrays[i], d1.Length);
            }


            for (int i1 = 0; i1< test_return_arrayofarrays.GetLength(0); i1++ )
            {
                double[] sub_Arr = test_return_arrayofarrays[i1];
                errorcode += "[";
                for (int i2 = 0; i2 < sub_Arr.GetLength(0); i2++)
                {
                    errorcode += test_return_arrayofarrays[i1][i2] + ",";

                }
                errorcode += "],";
            }
            //double[] d1 = (double [])arr.GetValue(0);
            //double[] d2 = (double[])arr.GetValue(1);
            //errorcode = "good " + d1[0] + " : " + d2[2];



            return true;
        }

        //[return: MarshalAs(UnmanagedType.SafeArray, SafeArraySubType = VarEnum.VT_VARIANT)]
        public double[,] TestGet2dArray()
        {
            double[,] arr;
            arr = new double[2, 2];
            arr[0, 0] = 1;
            return arr;
        }

        //[return: MarshalAs(UnmanagedType.SafeArray, SafeArraySubType = VarEnum.VT_VARIANT)]
        public string[] TestGetArrayOfArrays([MarshalAs(UnmanagedType.SafeArray, SafeArraySubType = VarEnum.VT_VARIANT)] Array arr)
        {
            int arr_num = arr.Length; //2 вложенных массива

            segment_list = new string[arr_num];

            for (int i = 0; i < arr_num; i++)  //по каждому массиву
            {
                double[] d1 = (double[])arr.GetValue(i); //получение i-массива
                segment_list[i] = GetDoubleArrayAsString(d1);
            }

            return segment_list;
        }
    }

}
