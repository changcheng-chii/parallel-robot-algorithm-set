using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Windows.Forms;
using Leadshine.SMC.IDE.Motion;
using System.Threading;
using System.Drawing.Drawing2D;
using MathNet.Numerics.LinearAlgebra;
using System.IO;

namespace 并联机床数控系统
{
    public class Inverse
    {
        public Inverse()
        { 

        }
        //尺寸参数极限范围
        double con;
        double[] a;   // 静平台半径;
        double[] b;   // 动平台半径;
        double d1min, d1max;  //连杆长度变化范围
        double d2min, d2max;  //连杆长度变化范围
        double d3min, d3max;  //连杆长度变化范围
        double R1_min, R1_max;//被动关节运动范围
        double R2_min, R2_max;
        double R3_min, R3_max;
        double U1_min, U1_max;   //静平台半径向量与支链夹角
        double U12_min, U12_max;
        double U2_min, U2_max;   //静平台半径向量与支链夹角
        double U22_min, U22_max;
        double S3_min, S3_max;   //静平台半径向量与支链夹角
        double e_tool;           //刀身长度，动平台原点到刀尖点距离，待确定
        
        //瞬时姿态记忆
        double psi, theta, z;                                           //并联动力头姿态参数
        public static double d1, d2, d3;                                         //支链长度
        double R1, R2, R3, U1, U2, U12, U22, S3;   //被动关节转角
        double[] vq;                                                 // 刀尖点现对于静系的坐标值
       
        //记录电机位置与机床位姿
        public static double[] PKM_X, PKM_Y, PKM_Z, PKM_A, PKM_B; //刀具路径相对于静坐标系坐标
        public static double[,] A1;                              //记录各位置下各运动轴位置
        public static double[,] B1;                             //记录所有当前位置到下一位置电机运动量
        public static double[,] D1;                            //记录所有当前位置到下一位置电机脉冲量
        public void inverse_run(double[] PKM_Pstar, double[] PKM_G, double[] G_X, double[] G_Y, double[] G_Z, double[] G_A, double[] G_B, double[] R)
        {
            //尺度确定
            con = Math.PI / 180;
            a = new double[3]; b = new double[3];
            a[0] = 系统界面.a1; a[1] = 系统界面.a2; a[2] = 系统界面.a3;   // 静平台半径;
            b[0] = 系统界面.b1; b[1] = 系统界面.b2; b[2] = 系统界面.b3;   // 动平台半径;
            d1min = 系统界面.d1min; d1max = 系统界面.d1max;  //连杆长度变化范围
            d2min = 系统界面.d2min; d2max = 系统界面.d2max;  //连杆长度变化范围
            d3min = 系统界面.d3min; d3max = 系统界面.d3max;  //连杆长度变化范围
            //关节运动范围
            R1_min = 系统界面.R1_min; R1_max = 系统界面.R1_max;
            R2_min = 系统界面.R2_min; R2_max = 系统界面.R2_max;
            R3_min = 系统界面.R3_min; R3_max = 系统界面.R3_max;

            U1_min = 系统界面.U1_min; U1_max = 系统界面.U1_max;     //静平台半径向量与支链夹角
            U12_min = 系统界面.U12_min; U12_max = 系统界面.U12_max;
            U2_min = 系统界面.U2_min; U2_max = 系统界面.U2_max;     //静平台半径向量与支链夹角
            U22_min = 系统界面.U22_min; U22_max = 系统界面.U22_max;
            S3_min = 系统界面.S3_min; S3_max = 系统界面.S3_max;     //静平台半径向量与支链夹角

            e_tool = 系统界面.e_tool;                                   //刀身长度，动平台原点到刀尖点距离，待确定

            motor_pulse(PKM_Pstar, PKM_G, G_X, G_Y, G_Z, G_A, G_B, R); 
        }
        double PKM_X0; double PKM_Y0; double PKM_Z0;
        double PKM_GX; double PKM_GY; double PKM_GZ;
        public void motor_pulse(double[] PKM_Pstar, double[] PKM_G, double[] G_X, double[] G_Y, double[] G_Z, double[] G_A, double[] G_B, double[] R)
        {
            PKM_X0 = PKM_Pstar[0]; PKM_Y0 = PKM_Pstar[1]; PKM_Z0 = PKM_Pstar[2];
            PKM_GX = PKM_G[0]; PKM_GY = PKM_G[1]; PKM_GZ = PKM_G[2];
            A1 = new double[G_Z.Length + 1, 5]; B1 = new double[G_Z.Length, 5]; D1 = new double[G_Z.Length, 6];  //G_z.length为D1有多少行？
            psi = 0;
            theta = 0;
            psi = psi * con; theta = theta * con;                   //角度转弧度
            z = PKM_Z0 - e_tool * Math.Cos(psi) * Math.Cos(theta);
            inverse_position(psi, theta, z);                        //调用运动学逆解
            if (R1_min < R1 && R1 < R1_max && R2_min < R2 && R2 < R2_max && R3_min < R3 && R3 < R3_max && U1_min < U1 && U1 < U1_max && U12_min < U12
                    && U12 < U12_max && U2_min < U2 && U2 < U2_max && U22_min < U22 && U22 < U22_max && S3_min < S3 && S3 < S3_max
                    && d1 < d1max && d2 < d2max && d3 < d3max && d1 > d1min && d2 > d2min && d3 > d3min) //边界条件
            {
                A1[0, 0] = d1;
                A1[0, 1] = d2;
                A1[0, 2] = d3;                        //并联动力头，并将电机初始位置设置为该位置           
                A1[0, 3] = -PKM_X0;
                A1[0, 4] = -PKM_Y0;                  //X-Y，并将电机初始位置设置为该位置

            }
            else
            {
                MessageBox.Show("所给刀尖初始点不在工作空间内，请核实初始点");
                return;
            }
            PKM_X = new double[G_Z.Length]; PKM_Y = new double[G_Z.Length]; PKM_Z = new double[G_Z.Length];
            PKM_A = new double[G_Z.Length]; PKM_B = new double[G_Z.Length];
            for (int i = 0; i < G_Z.Length; i = i + 1)   //将刀尖点轨迹从工件坐标系转到并联机架坐标系
            {
                PKM_X[i] = G_X[i]; PKM_Y[i] = G_Y[i]; PKM_Z[i] = PKM_GZ - G_Z[i];          //2019年04月30日 负号-   工件Z坐标-text的Z
                PKM_A[i] = G_A[i] * con; PKM_B[i] = G_B[i] * con;//角度转弧度
            }
            for (int j = 0; j < G_Z.Length; j = j + 1)
            {
                psi = PKM_A[j];
                theta = PKM_B[j];
                z = PKM_Z[j] - e_tool * Math.Cos(psi) * Math.Cos(theta);// 动平台理论中心点   
                inverse_position(psi, theta, z);                 //调用运动学逆解
                if (R1_min < R1 && R1 < R1_max && R2_min < R2 && R2 < R2_max && R3_min < R3 && R3 < R3_max && U1_min < U1 && U1 < U1_max && U12_min < U12
                    && U12 < U12_max && U2_min < U2 && U2 < U2_max && U22_min < U22 && U22 < U22_max && S3_min < S3 && S3 < S3_max
                    && d1 < d1max && d2 < d2max && d3 < d3max && d1 > d1min && d2 > d2min && d3 > d3min) //边界条件
                {
                    A1[j + 1, 0] = d1;
                    A1[j + 1, 1] = d2;
                    A1[j + 1, 2] = d3;             //并联动力头

                    A1[j + 1, 3] = vq[0] - PKM_X[j];
                    A1[j + 1, 4] = vq[1] - PKM_Y[j];  //X-Y

                    B1[j, 0] = A1[j + 1, 0] - A1[j, 0];//各轴相对位移量
                    B1[j, 1] = A1[j + 1, 1] - A1[j, 1];
                    B1[j, 2] = A1[j + 1, 2] - A1[j, 2];
                    B1[j, 3] = A1[j + 1, 3] - A1[j, 3];
                    B1[j, 4] = A1[j + 1, 4] - A1[j, 4];

                }
                else
                {
                    MessageBox.Show("刀尖轨迹点不在工作空间内");
                    return;
                }
            }
            int lp1 = 4;     //并联动力头丝杠导程
            int lp2 = 5;     //移动平台丝杠导程
            int fp = 10000;  //控制器分辨率   原先是10016
            for (int n = 0; n < G_Z.Length; n = n + 1)
            {
                D1[n, 0] = Math.Round(B1[n, 0] * fp / lp1);
                D1[n, 1] = Math.Round(B1[n, 1] * fp / lp1);
                D1[n, 2] = Math.Round(B1[n, 2] * fp / lp1);
                D1[n, 3] = Math.Round(B1[n, 3] * fp / lp2);//原来是lp12019年04月26日
                D1[n, 4] = Math.Round(B1[n, 4] * fp / lp2);
                D1[n, 5] = Math.Round(R[n] * fp / lp2);
            }
        }
             public void inverse_position(double psi, double theta, double z)
        {
            var JZ = Matrix<double>.Build; //初始化一个矩阵构建对象
            var XL = Vector<double>.Build;   //初始化一个向量的构建对象         
            double phi = 0;
            double x = z * Math.Tan(theta);
            double y = 0;                     //牵连运动

            double[] RotateA = { 1, 0, 0, 0, Math.Cos(psi), -Math.Sin(psi), 0, Math.Sin(psi), Math.Cos(psi) }; //绕X转矩阵
            double[] RotateB = { Math.Cos(theta), 0, Math.Sin(theta), 0, 1, 0, -Math.Sin(theta), 0, Math.Cos(theta) }; //绕Y轴转矩阵
            double[] RotateC = { Math.Cos(phi), -Math.Sin(phi), 0, Math.Sin(phi), Math.Cos(phi), 0, 0, 0, 1 };  //绕Z轴转矩阵

            var RZ = JZ.Dense(3, 3, RotateC).Transpose();  //将一维数组转化为3X3的多维数组,transpose转置
            var RY = JZ.Dense(3, 3, RotateB).Transpose(); 
            var RX = JZ.Dense(3, 3, RotateA).Transpose();
            var R = RZ * RY * RX;            //姿态变换矩阵

            a = new double[3]; b = new double[3];
            a[0] = 系统界面.a1; a[1] = 系统界面.a2; a[2] = 系统界面.a3;   // 静平台半径;
            b[0] = 系统界面.b1; b[1] = 系统界面.b2; b[2] = 系统界面.b3;   // 动平台半径;

            var vb10 = XL.Dense(3); vb10[0] = b[0] * 0;    vb10[1] = b[0] * (-1); vb10[2] = b[0] * 0; //vb10 = XL.Dense(3); 生成一个vb10的向量1X3
            var vb20 = XL.Dense(3); vb20[0] = b[1] * 0;    vb20[1] = b[1] * (1);  vb20[2] = b[1] * 0;
            var vb30 = XL.Dense(3); vb30[0] = b[2] * (1);  vb30[1] = b[2] * 0;    vb30[2] = b[2] * 0;//动平台各关节相对于动坐标系的坐标
            
            var vb1 = R * vb10;
            var vb2 = R * vb20;
            var vb3 = R * vb30;     //动坐标系原点与动平台各关节所成向量相对于静坐标系的坐标

            var va1 = XL.Dense(3); va1[0] = a[0] * 0;    va1[1] = a[0] * (-1); va1[2] = a[0] * 0;
            var va2 = XL.Dense(3); va2[0] = a[1] * 0;    va2[1] = a[1] * (1);  va2[2] = a[1] * 0;
            var va3 = XL.Dense(3); va3[0] = a[2] * (1);  va3[1] = a[2] * 0;    va3[2] = a[2] * 0;    //静平台各关节相对于静坐标系的坐标

            var vp = XL.Dense(3);  vp[0] = x;  vp[1] = y;  vp[2] = z;  //动坐标系原点相对于静坐标系原点的位置向量
            var vc0 = XL.Dense(3); vc0[0] = 0; vc0[1] = 1; vc0[2] = 0; //动平台求解U铰当前转角的参考向量
            var vc = R * vc0;
            var vd = XL.Dense(3); vd[0] = 1; vd[1] = 0; vd[2] = 0;     //静平台求解U铰当前转角的参考向量
            var vz0 = XL.Dense(3); vz0[0] = 0; vz0[1] = 0; vz0[2] = 1; //刀具方向向量
            var vz = R * vz0;
            var vq0 =vp + e_tool * vz;                                 //刀尖点相对于静坐标系坐标

            var qvw1 = vp + vb1 - va1;
            var qvw2 = vp + vb2 - va2;
            var qvw3 = vp + vb3 - va3;    //对应支链上两关节中心所成向量（静到动）

             d1 = Math.Sqrt(qvw1[0] * qvw1[0] + qvw1[1] * qvw1[1] + qvw1[2] * qvw1[2]);
             d2 = Math.Sqrt(qvw2[0] * qvw2[0] + qvw2[1] * qvw2[1] + qvw2[2] * qvw2[2]);
             d3 = Math.Sqrt(qvw3[0] * qvw3[0] + qvw3[1] * qvw3[1] + qvw3[2] * qvw3[2]);     //支链长度（对应支链两关节中心的距离）

            var vw1 = qvw1 / d1;
            var vw2 = qvw2 / d2;
            var vw3 = qvw3 / d3;      //对应支链的单位方向向量

            R1 = Math.Acos((vb1[0] * vw1[0] + vb1[1] * vw1[1] + vb1[2] * vw1[2]) / b[0]); R1 = R1 / con;
            U1 = Math.Acos((va1[0] * vw1[0] + va1[1] * vw1[1] + va1[2] * vw1[2]) / a[0]); U1 = U1 / con;
            U12 = Math.Acos(vd[0] * vw1[0] + vd[1] * vw1[1] + vd[2] * vw1[2]);            U12 = U12 / con;             //支链1被动关节当前姿态对应的转角
            R2 = Math.Acos((vb2[0] * vw2[0] + vb2[1] * vw2[1] + vb2[2] * vw2[2]) / b[1]); R2 = R2 / con;
            U2 = Math.Acos((va2[0] * vw2[0] + va2[1] * vw2[1] + va2[2] * vw2[2]) / a[1]); U2 = U2 / con;
            U22 = Math.Acos(vd[0] * vw2[0] + vd[1] * vw2[1] + vd[2] * vw2[2]);            U22 = U22 / con;             //支链2被动关节当前姿态对应的转角

            R3 = Math.Acos((va3[0] * vw3[0] + va3[1] * vw3[1] + va3[2] * vw3[2]) / a[2]); R3 = R3 / con;
            S3 = Math.Acos((vb3[0] * vw3[0] + vb3[1] * vw3[1] + vb3[2] * vw3[2]) / b[2]); S3 = S3 / con;
            //S32 = Math.Acos(vc[0] * vw3[0] + vc[1] * vw3[1] + vc[2] * vw3[2]);            S32 = S32 / con;             //支链3被动关节当前姿态对应的转角

             vq = new double[3];
             vq[0] = vq0[0]; vq[1] = vq0[1]; vq[2] = vq0[2];     //并联动力头刀尖坐标    
        }

    }
}
