using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Filter
{
    public class Accel
    {
        public TimeSpan m_Time;
        public int m_Pin = 0;
        public int m_vRef = 0;
        public int m_fx = 0;
        public int m_fy = 0;
        public int m_fz = 0;
        public double m_gx = 0;
        public double m_gy = 0;
        public double m_gz = 0;
        public double m_ax = 0;
        public double m_ay = 0;
        public double m_az = 0;

        public Accel(int ax, int ay, int az)
        {
            m_ax = ax;
            m_ay = ay;
            m_az = az;
        }

        public Accel(String data)
        {
            parse(data);
        }

        public int parse(String data)
        {
            int status = -1;
            String[] datas = data.Split(' ');
            if (datas.Length == 12)
            {
                m_Time = TimeSpan.Parse(datas[0]);
                m_Pin = Convert.ToInt32(datas[1]);
                m_vRef = Convert.ToInt32(datas[3]);
                m_fx = Convert.ToInt32(datas[5]);
                m_fy = Convert.ToInt32(datas[6]);
                m_fz = Convert.ToInt32(datas[7]);
                m_ax = Convert.ToInt32(datas[9]);
                m_ay = Convert.ToInt32(datas[10]);
                m_az = Convert.ToInt32(datas[11]);
            }
            return status;
        }
    }
}
