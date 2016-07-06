
#include <math.h> 
#include <cmath>
#include <iostream>
#include <Windows.h>
#include <tchar.h>
#include <string>
using namespace std;

//����opencv�е�ͷ�ļ�
#include <cv.h>   
#include <highgui.h> 


//ȫ�ֱ���
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
int width,height,channel;//Ӱ�񳤶ȣ���ȣ�ͨ����
double sigma = 0.4;//��׼��
double ratLow = 0.5;//����ֵ
double ratHigh = 0.79;//����ֵ

CvCapture *capture = cvCreateFileCapture("D:\\VideoData\\camera.avi");  //����ָ��capture��ʹ��ָ��ṹ��CVCapture���ҽṹ��CVCapture��������Ƶ��Ϣ

//��������
void RGBtoGray();

void Sobel();
void Roberts();
void Prewitt();
void Laplacian();
void Canny();//Canny����

int MyRoberts();
int MySobel();
int MyPrewitt();
int MyLaplacian();
int MyCanny();

void CreatGauss(double sigma, double **pdKernel, int *pnWidowSize);//����һά��˹����
void GaussianSmooth();//��ͼ���˹ƽ��
void Grad();//������,���ݶ�
void NonmaxSuppress();//�Ǽ���ֵ���Ʊ���
void Hysteresis();//���ú���Ѱ�ұ߽����
void EstimateThreshold(int *pThrHigh, int *pThrLow);//ͳ�Ƶ�ǰͼ��ֱ��ͼ���ж���ֵ
void TraceEdge(int y, int x, int nThrLow);//�����߽����
void noTraceEdge(int y,int x,int nLowThd);//�ǵݹ�

/*
   ������
   input:None
   output:return 0 
*/
int main()  
{  
	int enter;
	cout<<"Performing Edge Detection Algorithm"<<endl; //��ʾ����������л��˳�
	Sleep(2000);//��ͣ����
	system("cls");//����
	cout<<"Please Choose the Edge Detection Algorithm��"<<"\n"<<endl;
	cout<<"1��Roberts Edge Detection Algorithm "<<endl;//Roberts����
	cout<<"2��Prewitt Edge Detection Algorithm "<<endl;//Prewitt����
	cout<<"3��Sobel Edge Detection Algorithm "<<endl;//Sobel����
	cout<<"4��Laplacian Edge Detection Algorithm "<<endl;//Laplacian����
	cout<<"5��Canny Edge Detection Algorithm "<<endl;//Canny����
	
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
		frame = cvQueryFrame(capture);//���ļ���ץȡһ֡Ȼ���ѹ��������һ֡��
		if (!frame)
			break;
		width = (*frame).width;//���
		height =(*frame).height;//�߶�
		channel = (*frame).nChannels;//ͨ����
		m_pframe = new unsigned char [(*frame).height * (*frame).width * (*frame).nChannels];
		memcpy(m_pframe,(*frame).imageData,(*frame).imageSize);//����Դ�ڴ棨src��ָ����ڴ����򣩿�����Ŀ���ڴ�
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
	RGBtoGray();        /////�ҶȻ�
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
		frame = cvQueryFrame(capture);//���ļ���ץȡһ֡Ȼ���ѹ��������һ֡��
		if (!frame)
			break;
		width = (*frame).width;//���
		height =(*frame).height;//�߶�
		channel = (*frame).nChannels;//ͨ����
		m_pframe = new unsigned char [(*frame).height * (*frame).width * (*frame).nChannels];
		memcpy(m_pframe,(*frame).imageData,(*frame).imageSize);//����Դ�ڴ棨src��ָ����ڴ����򣩿�����Ŀ���ڴ�
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
	RGBtoGray();        //�ҶȻ�
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
		frame = cvQueryFrame(capture);//���ļ���ץȡһ֡Ȼ���ѹ��������һ֡��
		if (!frame)
			break;
		width = (*frame).width;//���
		height =(*frame).height;//�߶�
		channel = (*frame).nChannels;//ͨ����
		m_pframe = new unsigned char [(*frame).height * (*frame).width * (*frame).nChannels];
		memcpy(m_pframe,(*frame).imageData,(*frame).imageSize);//����Դ�ڴ棨src��ָ����ڴ����򣩿�����Ŀ���ڴ�
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
	RGBtoGray();//�ҶȻ�
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
		frame = cvQueryFrame(capture);//���ļ���ץȡһ֡Ȼ���ѹ��������һ֡��
		if (!frame)
			break;
		width = (*frame).width;//���
		height =(*frame).height;//�߶�
		channel = (*frame).nChannels;//ͨ����
		m_pframe = new unsigned char [(*frame).height * (*frame).width * (*frame).nChannels];
		memcpy(m_pframe,(*frame).imageData,(*frame).imageSize);//����Դ�ڴ棨src��ָ����ڴ����򣩿�����Ŀ���ڴ�
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
	RGBtoGray();//�ҶȻ�
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
	Canny����ģ��
	input:�Ҷ�ͼ��
	output����ģ����о�����ͼ��
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
	GaussianSmooth();//��˹�˲�
	Grad();//���ݶ�
	NonmaxSuppress();//�Ǽ���ֵ����
	Hysteresis();//˫��ֵ���
	delete []m_pGradX;//�ͷ�ָ��
	delete []m_pGradY;
	delete []m_pMag;
	for(int i = 0; i < height; i++)
		for(int j = 0; j < width; j++)
		{
			pDstFrame[i*width+j] =(char)m_pResult[i*width+j];//������
		}

}


//Canny����
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
		frame = cvQueryFrame(capture);//���ļ���ץȡһ֡Ȼ���ѹ��������һ֡��
		if (!frame)
			break;
		cvShowImage("Original_Video", frame);

		width = (*frame).width;//���
		height =(*frame).height;//�߶�
		channel = (*frame).nChannels;//ͨ����
		m_pframe = new unsigned char [(*frame).height * (*frame).width * (*frame).nChannels];
		memcpy(m_pframe,(*frame).imageData,(*frame).imageSize);//����Դ�ڴ棨src��ָ����ڴ����򣩿�����Ŀ���ڴ�

		Canny();//canny����

		IplImage* frame_Canny = cvCreateImage(cvGetSize(frame),IPL_DEPTH_8U,1);
		(*frame_Canny).imageData = pDstFrame;
		//cvCanny(Dst_iplFrame,cCanny,50,150,3.0);//�����Դ���cvCanny�㷨
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
	int nCenter;//�������ĵ�
	double dDis;//������һ�㵽���ĵ����
	//�м����
	double dValue;
	double dSum;
	long i;
	dSum = 0;
	//[-3*sigma,3*sigma] �������ݣ��Ḳ�Ǿ��󲿷��˲�ϵ��
	//dim = 1 + 2 * ((int) (3.0 * sigma))
	*pnWidowSize = 1+ 2*ceil(3*sigma);

	nCenter = (*pnWidowSize)/2;
	*pdKernel = new double[*pnWidowSize];

	//���ɸ�˹����
	for(i=0;i<(*pnWidowSize);i++)
	{
		dDis = double(i - nCenter);
		/*
		        1
		K=    ������  e ^ [-(x*x)/2��/��]
			  /2�Ц� 
		*/
		dValue = exp(-(0.5)*dDis*dDis/(sigma*sigma))/(sqrt(2*3.1415926)*sigma);
		(*pdKernel)[i] = dValue;
		dSum+=dValue;

	}
	//��һ��
	for(i=0;i<(*pnWidowSize);i++)
	{
		(*pdKernel)[i]/=dSum;
	}

}

void GaussianSmooth()//��˹�˲�
{
	long x,y,i;
	int windowSize;//��˹�˲�������
	int nlen; //���ڳ���

	double *pKernel; //һά��˹�˲���
	double dotMut; //��˹ϵ����ͼ�����ݵĵ��
	double weightSum; //�˲�ϵ���ܺ�
	double *pTemp;
	pTemp = new double[width*height];
	m_pfilter = new unsigned char[width*height];
	//����һά��˹����
	CreatGauss(sigma, &pKernel, &windowSize);//����һά��˹
	nlen = windowSize/2;

	//X�����˲�
	for(y = 0; y < height; y++)
		for(x= 0; x < width; x++)
		{
			//��ʼ������
			dotMut = 0;
			weightSum = 0;
			for(i = -nlen; i <= nlen; i++)
			{   //�ж��Ƿ���ͼ���ڲ�
				if((i+x)>=0 && (i+x)<width)
				{
					dotMut+=(double)m_pgray[y*width+(i+x)] * pKernel[nlen+i];
					weightSum += pKernel[nlen+i];
				}
			}

			pTemp[y*width+x] = dotMut/weightSum;
		}

	//Y�����˲�
	for(x= 0; x < width; x++)                      
		for(y = 0; y < height; y++)
		{
			//��ʼ������
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

// ������,���ݶ�
/*
	Sx = [ -1 1      Sy = [  1  1
	       -1 1 ]           -1 -1 ]
*/
void Grad()
{
	 long x,y;
	 //x����ķ�����
	 for(y=1;y<height-1;y++)
	 {
		 for(x=1;x<width-1;x++)
		 {
			 m_pGradX[y*width +x] = (int)( m_pfilter[y*width+x+1]-m_pfilter[y*width+ x-1]  );
		 }
	 }

	 //y��������
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
			 //���׷������ݶ�
			 dSqt1 = m_pGradX[y*width + x]*m_pGradX[y*width + x];//��X�����ݶ�
			 dSqt2 = m_pGradY[y*width + x]*m_pGradY[y*width + x];//��Y�����ݶȣ�
			 m_pMag[y*width+x] = (int)(sqrt(dSqt1+dSqt2)+0.5);//+0.5��������
		 }
	 }

}

//�Ǽ���ֵ���Ʊ���
void NonmaxSuppress()
{
	long y,x;
	int nPos;
	int gx,gy;//�ݶȷ���
	int g1,g2,g3,g4;//�м����
	double weight;
	double dTmp,dTmp1,dTmp2;

	//����ͼ���ԵΪ�����ܵķֽ��
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
			//��ǰ��
			nPos = y*width + x;
			//�����ǰ�����ݶȷ���Ϊ0�����Ǳ߽��
			if(m_pMag[nPos] == 0)
			{
				m_pResult[nPos] = 0;
			}
			else
			{
				//��ǰ����ݶȷ���
				dTmp = m_pMag[nPos];

				//x,y������
				gx = m_pGradX[nPos];
				gy = m_pGradY[nPos];

				//���������y������x������˵����������������y����
				if(abs(gy) > abs(gx))
				{
					//�����ֵ����
					weight = abs(gx)/abs(gy);

					g2 = m_pMag[nPos-width];
					g4 = m_pMag[nPos+width];

					//���x,y�����������ķ�����ͬ
                    //C Ϊ��ǰ���أ���g1-g4 ��λ�ù�ϵΪ��
                    //g1   g2
                    //     C
                    //     g4   g3
					if(gx*gy>0)
					{
						g1 = m_pMag[nPos-width-1];
						g3 = m_pMag[nPos+width+1];
					}
					//���x,y��������ķ����������෴
                    //C�ǵ�ǰ���أ���g1-g4�Ĺ�ϵΪ��
                    //      g2   g1
                    //      C
                    //g3   g4
					else
					{
						g1 = m_pMag[nPos-width+1];
						g3 = m_pMag[nPos+width-1];
					}

				}
				//���������x������y������˵�������ķ���������x����
				else
				{
					weight = abs(gy)/abs(gx);

					g2 = m_pMag[nPos+1];
					g4 = m_pMag[nPos-1];

					//���x,y��������ķ�����������ͬ
                    //��ǰ����C�� g1-g4�Ĺ�ϵΪ
                    //        g3
                    // g4  C  g2
                    // g1
					if(gx * gy > 0)
					{
						g1 = m_pMag[nPos+width+1];
						g3 = m_pMag[nPos-width-1];
					}
					//���x,y�����������ķ����෴
                    // C��g1-g4�Ĺ�ϵΪ
                    //        g1
                    // g4  C  g2
                    // g3
					{
						g1 = m_pMag[nPos-width+1];
						g3 = m_pMag[nPos+width-1];
					}

				}
				
				//weight = fabs(gx)/fabs(gy) = ctan(theta), thetaΪ�ݶȷ���.
				//weight = |dTemp1-g2|/|g1-g2| = |dTemp1-g2|/|C-g2| = ctan(theta); 

				dTmp1 = weight*g1 + (1-weight)*g2;
				dTmp2 = weight*g3 + (1-weight)*g4;

				//��ǰ���ص��ݶ��Ǿֲ������ֵ
				//�õ�����Ǳ߽��
				if(dTmp>=dTmp1 && dTmp>=dTmp2)
				{
					m_pResult[nPos] = 128;//�ǵ�ǰ��ĻҶ�ֵ
				}
				else
				{
					//�������Ǳ߽��
					m_pResult[nPos] = 0;
				}
			
			}
		}
	
}



//�ж���ֵ
void EstimateThreshold(int *pThrHigh, int *pThrLow)
{
	long y,x,k;
	int nHist[256]; //���ܱ߽���
	int nEdgeNum;//����ݶ���
	int nMaxMag = 0;
	int nHighCount;
	//��ʼ��
	for(k=0;k<256;k++)
	{
		nHist[k] = 0;
	}
	//ͳ��ֱ��ͼ,����ֱ��ͼ������ֵ
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

	//ͳ�ƾ����������ֵ���ơ����ж�������
	for(k=1;k<256;k++)
	{
		if(nHist[k] != 0)
		{
			nMaxMag = k;
		}

		//�ݶ�Ϊ0�ĵ��ǲ�����Ϊ�߽���
		//����non-maximum suppression���ж�������
		nEdgeNum += nHist[k];

	}
	//�ݶȱȸ���ֵ*pThrHigh С�����ص�����Ŀ
	nHighCount = (int)(ratHigh * nEdgeNum + 0.5);
	k=1;
	nEdgeNum = nHist[1];
	//�������ֵ
	while((k<(nMaxMag-1)) && (nEdgeNum < nHighCount))
	{
		k++;
		nEdgeNum += nHist[k];
	}
	*pThrHigh = k;
	//����ֵ=0.4����ֵ
	*pThrLow = (int)((*pThrHigh) * ratLow + 0.5);//��������

}

//�����߽����
void TraceEdge(int y, int x, int nThrLow)
{	 
	//   ��8�������ؽ��в�ѯ
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
		//   ���������Ϊ���ܵı߽��,��û�д����,�����ݶȴ�����ֵ 
		if(m_pResult[y_y * width + x_x] == 128 && m_pMag[y_y * width + x_x] >= nThrLow )
		{
			//�õ���Ϊ�߽��
			m_pResult[y_y * width + x_x] = 255;

			//�Ըõ�Ϊ�����ٽ��и���
			//TraceEdge(y_y,x_x,nThrLow);
		TraceEdge(y_y,x_x,nThrLow);    
		}
	}
}

//�ǵݹ����
void   noTraceEdge   (int  y,   int  x,   int  nLowThd)    
	{    
		//   ��8�������ؽ��в�ѯ  
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
				//   ���������Ϊ���ܵı߽�㣬��û�д����  
				//   �����ݶȴ�����ֵ  
				if(m_pResult[y_y*width+x_x]   ==   128   &&   m_pMag[y_y*width+x_x]>=nLowThd)  
				{  
					change=true;
					//   �Ѹõ����ó�Ϊ�߽��  
					m_pResult[y_y*width+x_x]   =   255   ;
					y=y_y;
					x=x_x;
					break;
					//   �Ըõ�Ϊ���Ľ��и���  
					//TraceEdge(yy,   xx,   nLowThd,   pUnchEdge,   pnMag,   nWidth);  
				}  
			}  
		}
	} 



//���ú���Ѱ�ұ߽����
void Hysteresis()
{
	long y,x;
	int nThrHigh,nThrLow;
	int nPos;

	EstimateThreshold(&nThrHigh,&nThrLow);//�ж���ֵ

	for(y=0;y<height;y++)
	{
		for(x=0;x<width;x++)
		{
			nPos = y*width + x;

			//����������ǿ��ܵı߽�㣬�����ݶȴ��ڸ���ֵ��
			//��������Ϊһ���߽�����
			if((m_pResult[nPos]==128) && (m_pMag[nPos] >= nThrHigh))
			{
				//���øõ�Ϊ�߽��
				m_pResult[nPos] = 255;
				noTraceEdge(y,x,nThrLow);
			}

		}
	}

	//�������Ѿ�������Ϊ�߽��
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

  



