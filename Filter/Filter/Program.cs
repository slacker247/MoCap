using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using System.Collections;

namespace Filter
{
    class Program
    {
        public static ArrayList m_Data;

        static void Main(String[] args)
        {
            m_Data = new ArrayList();

            // read three lines
            if (args.Length > 0)
            {
                System.Console.WriteLine("Loading Data...");
                loadData(args[0]);
                System.Console.WriteLine("...Loaded Data");
                // Convert accel data to gForce vectors
                System.Console.WriteLine("Calculating gForces...");
                calcGForces();
                System.Console.WriteLine("...Calculated gForces");
                // Smooth gForce vectors
                System.Console.WriteLine("Smoothing gForces...");
                smoothGForces();
                System.Console.WriteLine("...Smoothing gForces");
                // Convert gForce vectors to angles
                System.Console.WriteLine("gForces to Angles...");
                gForceToAngles();
                System.Console.WriteLine("...gForces to Angles");
                // Smooth angles
                System.Console.WriteLine("Smoothing Angles...");
                smoothAngles();
                System.Console.WriteLine("...Smoothing Angles");
                // Filter angles
                //averageForces(args[0]);
                writeData(args[0]);
                System.Console.WriteLine("Wrote Data.");
                //System.Console.ReadKey();
            }
        }

        public static void writeData(String fileName)
        {
            FileStream file = new FileStream(fileName + ".xaf", FileMode.OpenOrCreate, FileAccess.Write);
            StreamWriter sw = new StreamWriter(file);
            long start = ((Accel)m_Data[0]).m_Time.Ticks;
            int index = m_Data.Count - 1;

            bool test = true;
            while (test)
            {
                test = false;
                if(((Accel)m_Data[index]).m_Time.Ticks == 0)
                    test = true;
                if (double.IsNaN(((Accel)m_Data[index]).m_ax))
                    test = true;
                if (double.IsNaN(((Accel)m_Data[index]).m_ay))
                    test = true;
                if (double.IsNaN(((Accel)m_Data[index]).m_az))
                    test = true;
                index--;
            }

            TimeSpan sec = new TimeSpan(0, 0, 1);

            // TODO : Support multiple data  points - pin is the identifier.
            sw.WriteLine("<?xml version=\"1.0\" encoding=\"utf-8\"?>\n" +
                         "<MaxAnimation version=\"1.00\" date=\"Tue Jun 15 19:36:36 2010\">\n" +
                         "<SceneInfo fileName=\"\" startTick=\"0\" endTick=\"" + (((Accel)m_Data[index]).m_Time.Ticks - start) + "\" frameRate=\"30\" ticksPerFrame=\"" + (sec.Ticks/30) + "\" />\n" +
                         "<Node name=\"Sphere01\" parentNode=\"Scene Root\" parentNodeIndex=\"0\" numChildren=\"0\">\n" +
                         "<Controller name=\"Sphere01 \\ Transform \\ Rotation \\ Data \\ X Rotation\" classOf=\"Bezier Float\" classID=\"2007,0\" superClassID=\"9003\" subNum=\"0\" numChildren=\"0\" filterType=\"rotx\" outOfRangeBefore=\"constant\" outOfRangeAfter=\"constant\">\n" +
                         "<Keys count=\"" + (index + 1) + "\" inRangeLoop=\"false\" outRangeLoop=\"false\">");
            for (int i = 0; i <= index; i++)
            {
                Accel acl = ((Accel)m_Data[i]);
                sw.WriteLine("<Key t=\"" + (acl.m_Time.Ticks - start) + "\" xB=\"false\" yB=\"false\" zB=\"false\" wB=\"false\" cVel=\"false\" unconHan=\"false\" inTan=\"flat\" outTan=\"flat\" v=\"" + acl.m_ax + " \" inTanVal=\"0.0 \" outTanVal=\"0.0 \" inLen=\"0.0 \" outLen=\"0.0 \" />");
            }
            sw.WriteLine("</Keys>\n" +
                         "</Controller>\n" +
                         "<Controller name=\"Sphere01 \\ Transform \\ Rotation \\ Data \\ Y Rotation\" classOf=\"Bezier Float\" classID=\"2007,0\" superClassID=\"9003\" subNum=\"1\" numChildren=\"0\" filterType=\"roty\" outOfRangeBefore=\"constant\" outOfRangeAfter=\"constant\">\n" +
                         "<Keys count=\"" + (index + 1) + "\" inRangeLoop=\"false\" outRangeLoop=\"false\">");
            for (int i = 0; i <= index; i++)
            {
                Accel acl = ((Accel)m_Data[i]);
                sw.WriteLine("<Key t=\"" + (acl.m_Time.Ticks - start) + "\" xB=\"false\" yB=\"false\" zB=\"false\" wB=\"false\" cVel=\"false\" unconHan=\"false\" inTan=\"flat\" outTan=\"flat\" v=\"" + acl.m_ay + " \" inTanVal=\"0.0 \" outTanVal=\"0.0 \" inLen=\"0.0 \" outLen=\"0.0 \" />");
            }
            sw.WriteLine("</Keys>\n" +
                         "</Controller>\n" +
                         "<Controller name=\"Sphere01 \\ Transform \\ Rotation \\ Data \\ Z Rotation\" classOf=\"Bezier Float\" classID=\"2007,0\" superClassID=\"9003\" subNum=\"2\" numChildren=\"0\" filterType=\"rotz\" outOfRangeBefore=\"constant\" outOfRangeAfter=\"constant\">\n" +
                         "<Keys count=\"" + (index + 1) + "\" inRangeLoop=\"false\" outRangeLoop=\"false\">");
            for (int i = 0; i <= index; i++)
            {
                Accel acl = ((Accel)m_Data[i]);
                sw.WriteLine("<Key t=\"" + (acl.m_Time.Ticks - start) + "\" xB=\"false\" yB=\"false\" zB=\"false\" wB=\"false\" cVel=\"false\" unconHan=\"false\" inTan=\"flat\" outTan=\"flat\" v=\"" + acl.m_az + " \" inTanVal=\"0.0 \" outTanVal=\"0.0 \" inLen=\"0.0 \" outLen=\"0.0 \" />");
            }
            sw.WriteLine("</Keys>\n" +
                         "</Controller>\n" +
                         "</Node>\n" +
                         "</MaxAnimation>");

            sw.Close();
            file.Close();
        }

        public static void gForceToAngles()
        {
            for (int i = 0; i < m_Data.Count; i++)
            {
                Accel acl = ((Accel)m_Data[i]);
                acl.m_ax = angle(acl.m_gy, acl.m_gz);
                acl.m_ay = angle(acl.m_gx, acl.m_gz);
                acl.m_az = angle(acl.m_gx, acl.m_gy);

            }
         }

        public static double angle(double lhs, double rhs)
        {
            double ang = 0.0;
            double pDotQ = 0.0;
            double magPQ = 0.0;
            double aRads = 0.0;
    
            pDotQ = (lhs * 1) + (rhs * 0);
            magPQ = Math.Sqrt((Math.Pow(lhs, 2) + Math.Pow(rhs, 2)) * (Math.Pow(1, 2) + Math.Pow(0, 2)));
            aRads = Math.Acos(pDotQ / magPQ);
            ang = aRads;
            return ang;
        }

        public static void smoothAngles()
        {
            for (int i = 0; i < m_Data.Count; i++)
            {
                Accel acl = (Accel)m_Data[i];
                double sumx = 0;
                double sumy = 0;
                double sumz = 0;
                int x = 0;
                for (; x < 40; x++)
                {
                    if (i + x < m_Data.Count)
                    {
                        sumx += ((Accel)m_Data[i + x]).m_ax;
                        sumy += ((Accel)m_Data[i + x]).m_ay;
                        sumz += ((Accel)m_Data[i + x]).m_az;
                    }
                    else
                        break;
                }
                if (x != 0)
                {
                    acl.m_ax = sumx / (double)x;
                    acl.m_ay = sumy / (double)x;
                    acl.m_az = sumz / (double)x;
                }
            }
        }

        public static void smoothGForces()
        {
            for (int i = 0; i < m_Data.Count; i++)
            {
                Accel acl = (Accel)m_Data[i];
                double sumx = 0;
                double sumy = 0;
                double sumz = 0;
                int x = 0;
                for (; x < 40; x++)
                {
                    if (i + x < m_Data.Count)
                    {
                        sumx += ((Accel)m_Data[i + x]).m_gx;
                        sumy += ((Accel)m_Data[i + x]).m_gy;
                        sumz += ((Accel)m_Data[i + x]).m_gz;
                    }
                    else
                        break;
                }
                if (x != 0)
                {
                    acl.m_gx = sumx / (double)x;
                    acl.m_gy = sumy / (double)x;
                    acl.m_gz = sumz / (double)x;
                }
            }
        }

        public static void calcGForces()
        {
            for (int i = 0; i < m_Data.Count; i++)
            {
                Accel acl = (Accel)m_Data[i];
                acl.m_gx = (float)acl.m_fx / 455f;
                acl.m_gy = (float)acl.m_fy / 455f;
                acl.m_gz = (float)acl.m_fz / 455f;
            }
        }

        public static void filterForces(String fileName)
        {
        }

        public static void loadData(String fileName)
        {
            FileStream file = new FileStream(fileName, FileMode.Open, FileAccess.Read);
            StreamReader sr = new StreamReader(file);

            while (!sr.EndOfStream)
            {
                string s = sr.ReadLine();
                Accel a = new Accel(s);
                m_Data.Add(a);
            }

            sr.Close();
            file.Close();
        }

        public static void averageForces(String fileName)
        {
            float mil = (1f / 30f) * 1000f;
            TimeSpan window = new TimeSpan(0, 0, 0, 0, (int)mil);

            int index = 0;

            FileStream file = new FileStream(fileName + ".ft", FileMode.OpenOrCreate, FileAccess.Write);
            StreamWriter sw = new StreamWriter(file);

            while (index < m_Data.Count)
            {
                ArrayList temp = new ArrayList();
                TimeSpan start = ((Accel)m_Data[index]).m_Time;
                bool test = true;
                int offset = 1;
                while (test && (index + offset) < m_Data.Count)
                {
                    if ((((Accel)m_Data[index + offset]).m_Time - start) < window)
                    {
                        temp.Add(index + offset);
                    }
                    else
                        test = false;
                    offset++;
                }
                int sumX = 0;
                int sumY = 0;
                int sumZ = 0;
                int x = 0;

                for (; x < temp.Count; x++)
                {
                    //                        Accel t = ((Accel)m_Data[index]);
                    Accel t = (Accel)m_Data[(int)temp[x]];
                    sumX += t.m_fx;
                    sumY += t.m_fy;
                    sumX += t.m_fz;
                }
                //newData.Add(new Accel(sumX/x, sumY/x, sumZ/x));
                sw.WriteLine(sumX / x + " " + sumY / x + " " + sumZ / x);
                //sw.WriteLine(x2 + " " + y + " " + z);
                index += offset;
            }

            sw.Close();
            file.Close();
        }

        public static void filterAngle(String fileName)
        {
            //float mil = (1f / 30f) * 1000f;
            //TimeSpan window = new TimeSpan(0, 0, 0, 0, (int)mil);

            //int index = 0;

            //FileStream file = new FileStream(fileName + ".ft", FileMode.OpenOrCreate, FileAccess.Write);
            //StreamWriter sw = new StreamWriter(file);

            //while (index < m_Data.Count)
            //{
            //    ArrayList temp = new ArrayList();
            //    TimeSpan start = ((Accel)m_Data[index]).m_Time;
            //    bool test = true;
            //    int offset = 1;
            //    while (test && (index + offset) < m_Data.Count)
            //    {
            //        if ((((Accel)m_Data[index + offset]).m_Time - start) < window)
            //        {
            //            temp.Add(index + offset);
            //        }
            //        else
            //            test = false;
            //        offset++;
            //    }
            //    int sumX = 0;
            //    int sumY = 0;
            //    int sumZ = 0;
            //    int lastX = ((Accel)m_Data[0]).m_ax;
            //    int lastY = ((Accel)m_Data[0]).m_ay;
            //    int lastZ = ((Accel)m_Data[0]).m_az;
            //    int x = 0;

            //    for (; x < temp.Count; x++)
            //    {
            //        //                        Accel t = ((Accel)m_Data[index]);
            //        Accel t = (Accel)m_Data[(int)temp[x]];
            //        int x2 = 0, y = 0, z = 0;
            //        int clip = 290;
            //        if (Math.Abs(t.m_ax - lastX) < clip)
            //        {
            //            sumX += t.m_ax;
            //            x2 = t.m_ax;
            //        }
            //        else
            //        {
            //            sumX += t.m_ax - 360;
            //            x2 = t.m_ax - 360;
            //        }
            //        if (Math.Abs(t.m_ay - lastY) < clip)
            //        {
            //            sumY += t.m_ay;
            //            y = t.m_ay;
            //        }
            //        else
            //        {
            //            sumY += t.m_ay - 360;
            //            y = t.m_ay - 360;
            //        }
            //        if (Math.Abs(t.m_az - lastZ) < clip)
            //        {
            //            sumX += t.m_az;
            //            z = t.m_az;
            //        }
            //        else
            //        {
            //            sumZ += t.m_az - 360;
            //            z = t.m_az - 360;
            //        }
            //        lastX = x2;
            //        lastY = y;
            //        lastZ = z;
            //    }
            //    //newData.Add(new Accel(sumX/x, sumY/x, sumZ/x));
            //    sw.WriteLine(sumX/x + " " + sumY/x + " " + sumZ/x);
            //    //sw.WriteLine(x2 + " " + y + " " + z);
            //    index++;// += offset;
            //}

            //sw.Close();
            //file.Close();
        }
    }
}
