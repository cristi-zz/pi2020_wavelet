// OpenCVApplication.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "common.h"
#include<vector>

void testParcurgereSimplaDiblookStyle()
{
	char fname[MAX_PATH];
	while (openFileDlg(fname))
	{
		Mat src = imread(fname, CV_LOAD_IMAGE_GRAYSCALE);
		int height = src.rows;
		int width = src.cols;
		Mat dst = src.clone();

		double t = (double)getTickCount(); // Get the current time [s]

		// the fastest approach using the “diblook style”
		uchar *lpSrc = src.data;
		uchar *lpDst = dst.data;
		int w = (int) src.step; // no dword alignment is done !!!
		for (int i = 0; i<height; i++)
			for (int j = 0; j < width; j++) {
				uchar val = lpSrc[i*w + j];
				lpDst[i*w + j] = 255 - val;
			}

		// Get the current time again and compute the time difference [s]
		t = ((double)getTickCount() - t) / getTickFrequency();
		// Print (in the console window) the processing time in [ms] 
		printf("Time = %.3f [ms]\n", t * 1000);

		imshow("input image",src);
		imshow("negative image",dst);
		waitKey();
	}
}

std::vector<int> getLowVector(std::vector<int> fullVector)
{
	std::vector<int> v;
	int iStop = fullVector.size() / 2;
	for (int i = 0; i < iStop; i++)
	{
		int s = fullVector.at(2 * i) + fullVector.at(2 * i + 1);
		v.push_back(s / 2);
	}

	return v;
}

std::vector<int> getHighVector(std::vector<int> fullVector)
{
	std::vector<int> v;
	int iStop = fullVector.size() / 2;
	for (int i = 0; i < iStop; i++)
	{
		int s = fullVector.at(2 * i) - fullVector.at(2 * i + 1);
		v.push_back(s / 2);
	}

	return v;
}

std::vector<Mat_<uchar>> divideIntoFour(Mat_<uchar> originalImage)
{
	std::vector<Mat_<uchar>> results;
	int rows = originalImage.rows;
	int cols = originalImage.cols;

	Mat_<uchar> l = Mat_<uchar>(rows, cols / 2);
	Mat_<uchar> h = Mat_<uchar>(rows, cols / 2);

	for (int r = 0; r < rows; r++)
	{
		std::vector<int> half1Low;
		std::vector<int> half1High;
		for (int c = 0; c < cols; c++)
		{
			half1Low.push_back(originalImage(r, c));
			half1High.push_back(originalImage(r, c));
		}

		half1Low = getLowVector(half1Low);

		half1High = getHighVector(half1High);

		for (int c = 0; c < half1Low.size(); c++)
		{
			l(r, c) = half1Low.at(c);
			h(r, c) =  half1High.at(c);
		}
	}

	Mat_<uchar> ll = Mat_<uchar>(rows, cols);
	Mat_<uchar> lh = Mat_<uchar>(rows, cols);
	Mat_<uchar> hl = Mat_<uchar>(rows, cols);
	Mat_<uchar> hh = Mat_<uchar>(rows, cols);

	for (int c = 0; c < cols / 2; c++)
	{
		std::vector<int> half1Low;
		std::vector<int> half1High;
		for (int r = 0; r < rows; r++)
		{
			half1Low.push_back(l(r, c));
			half1High.push_back(l(r, c));
		}

		half1Low = getLowVector(half1Low);

		half1High = getHighVector(half1High);


		for (int r = 0; r < half1Low.size(); r++)
		{
			ll(r, c) = half1Low.at(r);
			lh(r, c) = half1High.at(r);
		}
	}

	for (int c = 0; c < cols / 2; c++)
	{
		std::vector<int> half1Low;
		std::vector<int> half1High;
		for (int r = 0; r < rows; r++)
		{
			half1Low.push_back(h(r, c));
			half1High.push_back(h(r, c));
		}

		half1Low = getLowVector(half1Low);

		half1High = getHighVector(half1High);

		for (int r = 0; r < half1Low.size(); r++)
		{
			hl(r, c) = half1Low.at(r);
			hh(r, c) = half1High.at(r);
		}
	}

	results.push_back(ll);
	results.push_back(lh);
	results.push_back(hl);
	results.push_back(hh);

	return results;
}

void testDecomposition()
{
	char fname[MAX_PATH];
	while (openFileDlg(fname))
	{
		Mat_<uchar> src = imread(fname, CV_LOAD_IMAGE_GRAYSCALE);
		std::cout << src.rows << ", " << src.cols << std::endl;
		std::vector<Mat_<uchar>> results = divideIntoFour(src);
		imshow("Original", src);
		imshow("LL", results.at(0));
		imshow("LH", results.at(1));
		imshow("HL", results.at(2));
		imshow("HH", results.at(3));

		waitKey(0);
	}
}

int main()
{
	testDecomposition();

	getchar();

	return 0;
}