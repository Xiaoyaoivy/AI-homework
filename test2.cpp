
#include <math.h> 
#include <cmath>
#include <iostream>
#include <Windows.h>
#include <tchar.h>
#include <string>
using namespace std;

//引入opencv中的头文件
#include <cv.h>   
#include <highgui.h> 


//全局变量
unsigned char *m_pframe;
unsigned char *m_pgray;
unsigned char *m_pResult;
unsigned char *m_pfilter;

//unsigned char *m_psobel;
//unsigned char *m_proberts;
//unsigned char *m_pprewitt;
//unsigned char *m_plaplacian;

char *pDstFrame,*Gray;
int *m_pGradX,*m_pGradY,*m_pMag;
int width,height,channel;//影像长度，宽度，通道数
double sigma = 0.4;//标准差
double ratLow = 0.5;//低阈值
double ratHigh = 0.79;//高阈值

CvCapture *capture = cvCreateFileCapture("D:\\VideoData\\camera.avi");  //创建指针capture，使其指向结构体CVCapture，且结构体CVCapture包含了视频信息

//函数声明
void RGBtoGray();

void Sobel();
void Roberts();
void Prewitt();
void Laplacian();
void Canny();//Canny算子

int MyRoberts();
int MySobel();
int MyPrewitt();
int MyLaplacian();
int MyCanny();

void CreatGauss(double sigma, double **pdKernel, int *pnWidowSize);//产生一维高斯数据
void GaussianSmooth();//对图像高斯平滑
void Grad();//方向导数,求梯度
void NonmaxSuppress();//非极大值抑制保留
void Hysteresis();//利用函数寻找边界起点
void EstimateThreshold(int *pThrHigh, int *pThrLow);//统计当前图像直方图，判定阈值
void TraceEdge(int y, int x, int nThrLow);//用来边界跟踪
void noTraceEdge(int y,int x,int nLowThd);//非递归

/*
   主函数
   input:None
   output:return 0 
*/
int main()  
{  
	int enter;
	cout<<"Performing Edge Detection Algorithm"<<endl; //提示程序继续运行或退出
	Sleep(2000);//暂停两秒
	system("cls");//清屏
	cout<<"Please Choose the Edge Detection Algorithm："<<"\n"<<endl;
	cout<<"1、Roberts Edge Detection Algorithm "<<endl;//Roberts算子
	cout<<"2、Prewitt Edge Detection Algorithm "<<endl;//Prewitt算子
	cout<<"3、Sobel Edge Detection Algorithm "<<endl;//Sobel算子
	cout<<"4、Laplacian Edge Detection Algorithm "<<endl;//Laplacian算子
	cout<<"5、Canny Edge Detection Algorithm "<<endl;//Canny算子
	
	cout<<"Enter your choice:";
	cin>>enter;
	while(enter!=1&&enter!=2&&enter!=3&&enter!=4&&enter!=5)
	 {
		 cin.clear();
		 cin.sync();
		 cout << "Not number or wrong number"<<endl;
		 cout << "Enter your choice:";
		 cin>>enter;
	 }

	switch (enter)
	{
		case 1:
		{
			cout<<"Roberts processing"<<endl;
			MyRoberts();break;
		}
		case 2:
		{
			cout<<"Prewitt processing"<<endl;
			MyPrewitt();break;
		}
		case 3:
		{
			cout<<"Sobel processing"<<endl;
			MySobel();break;
		}
		case 4:
		{
			cout<<"Laplacian processing"<<endl;
			MyLaplacian();break;
		}
		case 5:
		{
			cout<<"Canny processing"<<endl;
			MyCanny();break;
		}
	}
	
	  return 0;
}


void RGBtoGray()
{
	m_pgray = new unsigned char[height*width];
	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++)
		{
			double r = m_pframe[(i*width+j)*channel+2]; 
			double g = m_pframe[(i*width+j)*channel+1];
			double b = m_pframe[(i*width+j)*channel+0];
			double gray = (r*30 + g*59 + b*11 + 50) / 100;
			m_pgray[(i*width+j)] = gray;
		}
	}
}



int MySobel()
{
	cvNamedWindow("Original_Video", CV_WINDOW_AUTOSIZE);
	cvNamedWindow("Sobel_Video", CV_WINDOW_AUTOSIZE);
	if (capture == NULL)
    {
		return 0;
	}
	IplImage *frame = NULL;	
	while(1)
	{
		frame = cvQueryFrame(capture);//从文件中抓取一帧然后解压并返回这一帧。
		if (!frame)
			break;
		width = (*frame).width;//宽度
		height =(*frame).height;//高度
		channel = (*frame).nChannels;//通道数
		m_pframe = new unsigned char [(*frame).height * (*frame).width * (*frame).nChannels];
		memcpy(m_pframe,(*frame).imageData,(*frame).imageSize);//把资源内存（src所指向的内存区域）拷贝到目标内存
		Sobel();
		IplImage* frame_Sobel = cvCreateImage(cvGetSize(frame),IPL_DEPTH_8U,1);
		(*frame_Sobel).imageData = pDstFrame;
		cvShowImage("Original_Video", frame);
		cvShowImage("Sobel_Video",frame_Sobel);
		char c = cvWaitKey(33);
		if (c == 27)
		{
			break;
		}		
	}
	cvReleaseCapture(&capture);
	cvDestroyAllWindows();  	
	return 0;
}

void Sobel()
{
	pDstFrame = new char[width*height];
	RGBtoGray();        /////灰度化
	for(int i = 1; i < height-1; i++)
		for(int j = 1; j < width-1; j++)
		{
			char ux=-(char)m_pgray[(i-1)*width+j-1]-2*(char)m_pgray[(i-1)*width+j]-(char)m_pgray[(i-1)*width+j+1]+(char)m_pgray[(i+1)*width+j-1]+2*(char)m_pgray[(i+1)*width+j]+(char)m_pgray[(i+1)*width+j+1];
			char uy=-(char)m_pgray[(i-1)*width+j-1]+(char)m_pgray[(i-1)*width+j+1]-2*(char)m_pgray[i*width+j-1]+2*(char)m_pgray[i*width+j+1]-(char)m_pgray[(i+1)*width+j-1]+(char)m_pgray[(i+1)*width+j+1];
			pDstFrame[i*width+j] = sqrt(pow(ux,2.0) + pow(uy,2.0));
		}
}

int MyRoberts()
{
	cvNamedWindow("Original_Video", CV_WINDOW_AUTOSIZE);
	cvNamedWindow("Roberts_Video", CV_WINDOW_AUTOSIZE);
	if (capture == NULL)
    {
		return 0;
	}
	IplImage *frame = NULL;	
	while(1)
	{
		frame = cvQueryFrame(capture);//从文件中抓取一帧然后解压并返回这一帧。
		if (!frame)
			break;
		width = (*frame).width;//宽度
		height =(*frame).height;//高度
		channel = (*frame).nChannels;//通道数
		m_pframe = new unsigned char [(*frame).height * (*frame).width * (*frame).nChannels];
		memcpy(m_pframe,(*frame).imageData,(*frame).imageSize);//把资源内存（src所指向的内存区域）拷贝到目标内存
		Roberts();
		IplImage* frame_Roberts = cvCreateImage(cvGetSize(frame),IPL_DEPTH_8U,1);
		(*frame_Roberts).imageData = pDstFrame;
		cvShowImage("Original_Video", frame);
		cvShowImage("Roberts_Video",frame_Roberts);
		char c = cvWaitKey(33);
		if (c == 27)
		{
			break;
		}		
	}
	cvReleaseCapture(&capture);
	cvDestroyAllWindows();  	
	return 0;
}

void Roberts()
{
	pDstFrame = new char[width*height];
	RGBtoGray();        //灰度化
	for(int i = 1; i < height-1; i++)
		for(int j = 1; j < width-1; j++)
		{
			char u=abs((char)m_pgray[(i-1)*width+j-1]-(char)m_pgray[(i+1)*width+j+1])+abs((char)m_pgray[(i+1)*width+j]-(char)m_pgray[(i)*width+j+1]);
			pDstFrame[i*width+j] = u;
		}
}

int MyPrewitt()
{
	cvNamedWindow("Original_Video", CV_WINDOW_AUTOSIZE);
	cvNamedWindow("Prewitt_Video", CV_WINDOW_AUTOSIZE);
	if (capture == NULL)
    {
		return 0;
	}
	IplImage *frame = NULL;	
	while(1)
	{
		frame = cvQueryFrame(capture);//从文件中抓取一帧然后解压并返回这一帧。
		if (!frame)
			break;
		width = (*frame).width;//宽度
		height =(*frame).height;//高度
		channel = (*frame).nChannels;//通道数
		m_pframe = new unsigned char [(*frame).height * (*frame).width * (*frame).nChannels];
		memcpy(m_pframe,(*frame).imageData,(*frame).imageSize);//把资源内存（src所指向的内存区域）拷贝到目标内存
		Prewitt();
		IplImage* frame_Prewitt = cvCreateImage(cvGetSize(frame),IPL_DEPTH_8U,1);
		(*frame_Prewitt).imageData = pDstFrame;
		cvShowImage("Original_Video", frame);
		cvShowImage("Prewitt_Video",frame_Prewitt);
		char c = cvWaitKey(33);
		if (c == 27)
		{
			break;
		}		
	}
	cvReleaseCapture(&capture);
	cvDestroyAllWindows();  	
	return 0;
}

void Prewitt()
{
	pDstFrame = new char[width*height];
	RGBtoGray();//灰度化
	for(int i = 1; i < height-1; i++)
		for(int j = 1; j < width-1; j++)
		{
			char a00=(char)m_pgray[(i-1)*width+j-1];
			char a01=(char)m_pgray[(i-1)*width+j];
			char a02=(char)m_pgray[(i-1)*width+j+1];

			char a10=(char)m_pgray[i*width+j-1];
			char a11=(char)m_pgray[i*width+j];
			char a12=(char)m_pgray[i*width+j+1];
			
			char a20=(char)m_pgray[(i+1)*width+j-1];
			char a21=(char)m_pgray[(i+1)*width+j-1];
			char a22=(char)m_pgray[(i+1)*width+j-1];

			char ux = a00 * (-1) + a01 * (-1) + a02 * (-1) + a20 * (1) + a21 * (1) + a22 * (1);
			char uy = (a00 * (-1) + a10 * (-1) + a20 * (-1)) + a02 * (1) + a12 * (1) + a22 * (1);
			pDstFrame[i*width+j] = sqrt(pow(ux,2.0) + pow(uy,2.0));
		}
}

int MyLaplacian()
{
	cvNamedWindow("Original_Video", CV_WINDOW_AUTOSIZE);
	cvNamedWindow("Laplacian_Video", CV_WINDOW_AUTOSIZE);
	if (capture == NULL)
    {
		return 0;
	}
	IplImage *frame = NULL;	
	while(1)
	{
		frame = cvQueryFrame(capture);//从文件中抓取一帧然后解压并返回这一帧。
		if (!frame)
			break;
		width = (*frame).width;//宽度
		height =(*frame).height;//高度
		channel = (*frame).nChannels;//通道数
		m_pframe = new unsigned char [(*frame).height * (*frame).width * (*frame).nChannels];
		memcpy(m_pframe,(*frame).imageData,(*frame).imageSize);//把资源内存（src所指向的内存区域）拷贝到目标内存
		Laplacian();
		IplImage* frame_Laplacian = cvCreateImage(cvGetSize(frame),IPL_DEPTH_8U,1);
		(*frame_Laplacian).imageData = pDstFrame;
		cvShowImage("Original_Video", frame);
		cvShowImage("Laplacian_Video",frame_Laplacian);
		char c = cvWaitKey(33);
		if (c == 27)
		{
			break;
		}		
	}
	cvReleaseCapture(&capture);
	cvDestroyAllWindows();  	
	return 0;
}


void Laplacian()
{
	pDstFrame = new char[width*height];
	RGBtoGray();//灰度化
	for(int i = 1; i < height-1; i++)
		for(int j = 1; j < width-1; j++)
		{
			char a00=(char)m_pgray[(i-1)*width+j-1];
			char a01=(char)m_pgray[(i-1)*width+j];
			char a02=(char)m_pgray[(i-1)*width+j+1];

			char a10=(char)m_pgray[i*width+j-1];
			char a11=(char)m_pgray[i*width+j];
			char a12=(char)m_pgray[i*width+j+1];
			
			char a20=(char)m_pgray[(i+1)*width+j-1];
			char a21=(char)m_pgray[(i+1)*width+j-1];
			char a22=(char)m_pgray[(i+1)*width+j-1];

			double u = -1 * (a00 + a01 + a02 + a10 + a12 + a20 + a21 + a22) + 8 * a11;

			pDstFrame[i*width+j] = u;
		}
}


 /*
	Canny算子模版
	input:灰度图像
	output：与模版进行卷积后的图像
*/
void Canny()
{
	m_pGradX = new int[width*height];
	m_pGradY = new int[width*height];
	m_pMag = new int[width*height];
	pDstFrame = new char[width*height];
	m_pResult = new unsigned char[width*height];
	m_pgray = new unsigned char[height*width];
	RGBtoGray();
	GaussianSmooth();//高斯滤波
	Grad();//求梯度
	NonmaxSuppress();//非极大值抑制
	Hysteresis();//双阈值检测
	delete []m_pGradX;//释放指针
	delete []m_pGradY;
	delete []m_pMag;
	for(int i = 0; i < height; i++)
		for(int j = 0; j < width; j++)
		{
			pDstFrame[i*width+j] =(char)m_pResult[i*width+j];//输出结果
		}

}


//Canny算子
int MyCanny()
{
	cvNamedWindow("Original_Video", CV_WINDOW_AUTOSIZE);
	cvNamedWindow("Canny_Video", CV_WINDOW_AUTOSIZE);
	
	if (capture == NULL)
    {
		return 0;
	}

	IplImage *frame = NULL;
	int iNum=0, iFrameH=0, iFrameW=0, iFps=0, iNumFrames=0;

	char ch[10];
	string strImageName;

	while(1)
	{
		frame = cvQueryFrame(capture);//从文件中抓取一帧然后解压并返回这一帧。
		if (!frame)
			break;
		cvShowImage("Original_Video", frame);

		width = (*frame).width;//宽度
		height =(*frame).height;//高度
		channel = (*frame).nChannels;//通道数
		m_pframe = new unsigned char [(*frame).height * (*frame).width * (*frame).nChannels];
		memcpy(m_pframe,(*frame).imageData,(*frame).imageSize);//把资源内存（src所指向的内存区域）拷贝到目标内存

		Canny();//canny算子

		IplImage* frame_Canny = cvCreateImage(cvGetSize(frame),IPL_DEPTH_8U,1);
		(*frame_Canny).imageData = pDstFrame;
		//cvCanny(Dst_iplFrame,cCanny,50,150,3.0);//程序自带的cvCanny算法
		cvShowImage("Canny_Video",frame_Canny);

		char c = cvWaitKey(33);
		if (c == 27)
		{
			break;
		}
			
	}
	cvReleaseCapture(&capture);

	cvDestroyAllWindows();  
	
	return 0;
	
}

void CreatGauss(double sigma, double **pdKernel, int *pnWidowSize)
{
	int nCenter;//数组中心点
	double dDis;//数组中一点到中心点距离
	//中间变量
	double dValue;
	double dSum;
	long i;
	dSum = 0;
	//[-3*sigma,3*sigma] 以内数据，会覆盖绝大部分滤波系数
	//dim = 1 + 2 * ((int) (3.0 * sigma))
	*pnWidowSize = 1+ 2*ceil(3*sigma);

	nCenter = (*pnWidowSize)/2;
	*pdKernel = new double[*pnWidowSize];

	//生成高斯数据
	for(i=0;i<(*pnWidowSize);i++)
	{
		dDis = double(i - nCenter);
		/*
		        1
		K=    ———  e ^ [-(x*x)/2δ/δ]
			  /2πδ 
		*/
		dValue = exp(-(0.5)*dDis*dDis/(sigma*sigma))/(sqrt(2*3.1415926)*sigma);
		(*pdKernel)[i] = dValue;
		dSum+=dValue;

	}
	//归一化
	for(i=0;i<(*pnWidowSize);i++)
	{
		(*pdKernel)[i]/=dSum;
	}

}

void GaussianSmooth()//高斯滤波
{
	long x,y,i;
	int windowSize;//高斯滤波器长度
	int nlen; //窗口长度

	double *pKernel; //一维高斯滤波器
	double dotMut; //高斯系数与图像数据的点乘
	double weightSum; //滤波系数总和
	double *pTemp;
	pTemp = new double[width*height];
	m_pfilter = new unsigned char[width*height];
	//产生一维高斯数据
	CreatGauss(sigma, &pKernel, &windowSize);//创建一维高斯
	nlen = windowSize/2;

	//X方向滤波
	for(y = 0; y < height; y++)
		for(x= 0; x < width; x++)
		{
			//初始化参数
			dotMut = 0;
			weightSum = 0;
			for(i = -nlen; i <= nlen; i++)
			{   //判断是否在图像内部
				if((i+x)>=0 && (i+x)<width)
				{
					dotMut+=(double)m_pgray[y*width+(i+x)] * pKernel[nlen+i];
					weightSum += pKernel[nlen+i];
				}
			}

			pTemp[y*width+x] = dotMut/weightSum;
		}

	//Y方向滤波
	for(x= 0; x < width; x++)                      
		for(y = 0; y < height; y++)
		{
			//初始化参数
			dotMut = 0;
			weightSum = 0;
			for(i = -nlen; i <= nlen; i++)
			{
				if((i+y)>=0 && (i+y)<height)
				{
					dotMut+=(double)pTemp[(y+i)*width+x] * pKernel[nlen+i];
					weightSum += pKernel[nlen+i];
				}
			}

			m_pfilter[y*width+x] = (unsigned char)(int)dotMut/weightSum;
		}

	delete []pKernel;
	pKernel = NULL;
	delete []pTemp;
	pTemp = NULL;

}

// 方向导数,求梯度
/*
	Sx = [ -1 1      Sy = [  1  1
	       -1 1 ]           -1 -1 ]
*/
void Grad()
{
	 long x,y;
	 //x方向的方向导数
	 for(y=1;y<height-1;y++)
	 {
		 for(x=1;x<width-1;x++)
		 {
			 m_pGradX[y*width +x] = (int)( m_pfilter[y*width+x+1]-m_pfilter[y*width+ x-1]  );
		 }
	 }

	 //y方向方向导数
	 for(x=1;x<width-1;x++)
	 {
		 for(y=1;y<height-1;y++)
		 {
			 m_pGradY[y*width +x] = (int)(m_pfilter[(y+1)*width +x] - m_pfilter[(y-1)*width +x]);
		 }
	 }

	 double dSqt1;
	 double dSqt2;

	 for(y=0; y<height; y++)
	 {
		 for(x=0; x<width; x++)
		 {
			 //二阶范数求梯度
			 dSqt1 = m_pGradX[y*width + x]*m_pGradX[y*width + x];//求X方向梯度
			 dSqt2 = m_pGradY[y*width + x]*m_pGradY[y*width + x];//求Y方向梯度；
			 m_pMag[y*width+x] = (int)(sqrt(dSqt1+dSqt2)+0.5);//+0.5四舍五入
		 }
	 }

}

//非极大值抑制保留
void NonmaxSuppress()
{
	long y,x;
	int nPos;
	int gx,gy;//梯度分量
	int g1,g2,g3,g4;//中间变量
	double weight;
	double dTmp,dTmp1,dTmp2;

	//设置图像边缘为不可能的分界点
	for(x=0;x<width;x++)
	{
		m_pResult[x] = 0;
		m_pResult[(height-1)*width+x] = 0;

	}
	for(y=0;y<height;y++)
	{
		m_pResult[y*width] = 0;
		m_pResult[y*width + width-1] = 0;
	}

	for (y = 1; y < height-1; y++)
		for(x =1; x < width-1; x++)
		{
			//当前点
			nPos = y*width + x;
			//如果当前像素梯度幅度为0，则不是边界点
			if(m_pMag[nPos] == 0)
			{
				m_pResult[nPos] = 0;
			}
			else
			{
				//当前点的梯度幅度
				dTmp = m_pMag[nPos];

				//x,y方向导数
				gx = m_pGradX[nPos];
				gy = m_pGradY[nPos];

				//如果方向导数y分量比x分量大，说明导数方向趋向于y分量
				if(abs(gy) > abs(gx))
				{
					//计算插值比例
					weight = abs(gx)/abs(gy);

					g2 = m_pMag[nPos-width];
					g4 = m_pMag[nPos+width];

					//如果x,y两个方向导数的符号相同
                    //C 为当前像素，与g1-g4 的位置关系为：
                    //g1   g2
                    //     C
                    //     g4   g3
					if(gx*gy>0)
					{
						g1 = m_pMag[nPos-width-1];
						g3 = m_pMag[nPos+width+1];
					}
					//如果x,y两个方向的方向导数方向相反
                    //C是当前像素，与g1-g4的关系为：
                    //      g2   g1
                    //      C
                    //g3   g4
					else
					{
						g1 = m_pMag[nPos-width+1];
						g3 = m_pMag[nPos+width-1];
					}

				}
				//如果方向导数x分量比y分量大，说明导数的方向趋向于x分量
				else
				{
					weight = abs(gy)/abs(gx);

					g2 = m_pMag[nPos+1];
					g4 = m_pMag[nPos-1];

					//如果x,y两个方向的方向导数符号相同
                    //当前像素C与 g1-g4的关系为
                    //        g3
                    // g4  C  g2
                    // g1
					if(gx * gy > 0)
					{
						g1 = m_pMag[nPos+width+1];
						g3 = m_pMag[nPos-width-1];
					}
					//如果x,y两个方向导数的方向相反
                    // C与g1-g4的关系为
                    //        g1
                    // g4  C  g2
                    // g3
					{
						g1 = m_pMag[nPos-width+1];
						g3 = m_pMag[nPos+width-1];
					}

				}
				
				//weight = fabs(gx)/fabs(gy) = ctan(theta), theta为梯度方向.
				//weight = |dTemp1-g2|/|g1-g2| = |dTemp1-g2|/|C-g2| = ctan(theta); 

				dTmp1 = weight*g1 + (1-weight)*g2;
				dTmp2 = weight*g3 + (1-weight)*g4;

				//当前像素的梯度是局部的最大值
				//该点可能是边界点
				if(dTmp>=dTmp1 && dTmp>=dTmp2)
				{
					m_pResult[nPos] = 128;//是当前点的灰度值
				}
				else
				{
					//不可能是边界点
					m_pResult[nPos] = 0;
				}
			
			}
		}
	
}



//判定阈值
void EstimateThreshold(int *pThrHigh, int *pThrLow)
{
	long y,x,k;
	int nHist[256]; //可能边界数
	int nEdgeNum;//最大梯度数
	int nMaxMag = 0;
	int nHighCount;
	//初始化
	for(k=0;k<256;k++)
	{
		nHist[k] = 0;
	}
	//统计直方图,利用直方图计算阈值
	for(y=0;y<height;y++)
	{
		for(x=0;x<width;x++)
		{
			if(m_pgray[y*width+x]==128)
			{
				nHist[m_pMag[y*width+x]]++;
			}
		}
	}
	nEdgeNum = nHist[0];

	//统计经过“非最大值抑制”后有多少像素
	for(k=1;k<256;k++)
	{
		if(nHist[k] != 0)
		{
			nMaxMag = k;
		}

		//梯度为0的点是不可能为边界点的
		//经过non-maximum suppression后有多少像素
		nEdgeNum += nHist[k];

	}
	//梯度比高阈值*pThrHigh 小的像素点总数目
	nHighCount = (int)(ratHigh * nEdgeNum + 0.5);
	k=1;
	nEdgeNum = nHist[1];
	//计算高阈值
	while((k<(nMaxMag-1)) && (nEdgeNum < nHighCount))
	{
		k++;
		nEdgeNum += nHist[k];
	}
	*pThrHigh = k;
	//低阈值=0.4高阈值
	*pThrLow = (int)((*pThrHigh) * ratLow + 0.5);//四舍五入

}

//用来边界跟踪
void TraceEdge(int y, int x, int nThrLow)
{	 
	//   对8邻域象素进行查询
	int xNum[8] = { 1,  1,  0,
		           -1,     -1,
				   -1,  0,  1};
	int yNum[8] = { 0,  1,  1,
		            1,      0,
				   -1, -1, -1};
	long y_y, x_x, k;

	for(k = 0;k < 8;k ++)
	{
		y_y = y + yNum[k];
		x_x = x + xNum[k];
		//   如果该象素为可能的边界点,又没有处理过,并且梯度大于阈值 
		if(m_pResult[y_y * width + x_x] == 128 && m_pMag[y_y * width + x_x] >= nThrLow )
		{
			//该点设为边界点
			m_pResult[y_y * width + x_x] = 255;

			//以该点为中心再进行跟踪
			//TraceEdge(y_y,x_x,nThrLow);
		TraceEdge(y_y,x_x,nThrLow);    
		}
	}
}

//非递归调用
void   noTraceEdge   (int  y,   int  x,   int  nLowThd)    
	{    
		//   对8邻域象素进行查询  
		int   xNum[8]   =   {1,   1,   0,-1,-1,-1,   0,   1}   ;  
		int   yNum[8]   =   {0,   1,   1,   1,0   ,-1,-1,-1}   ;  

		int   y_y   ;  
		int   x_x   ;  
		int   k   ;  

		bool change=true;
		while(change)
		{
			change=false;
			for(k=0;   k<8;   k++)  
			{  
				y_y   =   y   +   yNum[k]   ;  
				x_x   =   x   +   xNum[k]   ;  
				//   如果该象素为可能的边界点，又没有处理过  
				//   并且梯度大于阈值  
				if(m_pResult[y_y*width+x_x]   ==   128   &&   m_pMag[y_y*width+x_x]>=nLowThd)  
				{  
					change=true;
					//   把该点设置成为边界点  
					m_pResult[y_y*width+x_x]   =   255   ;
					y=y_y;
					x=x_x;
					break;
					//   以该点为中心进行跟踪  
					//TraceEdge(yy,   xx,   nLowThd,   pUnchEdge,   pnMag,   nWidth);  
				}  
			}  
		}
	} 



//利用函数寻找边界起点
void Hysteresis()
{
	long y,x;
	int nThrHigh,nThrLow;
	int nPos;

	EstimateThreshold(&nThrHigh,&nThrLow);//判断阈值

	for(y=0;y<height;y++)
	{
		for(x=0;x<width;x++)
		{
			nPos = y*width + x;

			//如果该像素是可能的边界点，并且梯度大于高阈值，
			//该像素作为一个边界的起点
			if((m_pResult[nPos]==128) && (m_pMag[nPos] >= nThrHigh))
			{
				//设置该点为边界点
				m_pResult[nPos] = 255;
				noTraceEdge(y,x,nThrLow);
			}

		}
	}

	//其他点已经不可能为边界点
	for(y=0;y<height;y++)
	{
		for(x=0;x<width;x++)
		{
			nPos = y*width + x;

			if(m_pResult[nPos] != 255)
			{
				m_pResult[nPos] = 0;
			}
		}
	}

}

  



