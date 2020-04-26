// OpenCVApplication.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "common.h"
#include<vector>

int hVec[2] = { 1, -1 };
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

										   // the fastest approach using the �diblook style�
		uchar *lpSrc = src.data;
		uchar *lpDst = dst.data;
		int w = (int)src.step; // no dword alignment is done !!!
		for (int i = 0; i < height; i++)
			for (int j = 0; j < width; j++) {
				uchar val = lpSrc[i*w + j];
				lpDst[i*w + j] = 255 - val;
			}

		// Get the current time again and compute the time difference [s]
		t = ((double)getTickCount() - t) / getTickFrequency();
		// Print (in the console window) the processing time in [ms] 
		printf("Time = %.3f [ms]\n", t * 1000);

		imshow("input image", src);
		imshow("negative image", dst);
		waitKey();
	}
}

float clampValue(float val, float min, float max)
{
	if (val < min)
	{
		return min;
	}
	if (val > max)
	{
		return max;
	}
	return val;
}

std::vector<float> getLowVector(std::vector<float> fullVector)
{
	std::vector<float> v;
	int iStop = fullVector.size() / 2;
	for (int i = 0; i < iStop; i++)
	{
		float s = fullVector.at(2 * i) + fullVector.at(2 * i + 1);
		v.push_back(s / 2);
	}

	return v;
}

std::vector<float> getHighVector(std::vector<float> fullVector)
{
	std::vector<float> v;
	int iStop = fullVector.size() / 2;
	for (int i = 0; i < iStop; i++)
	{
		float s = fullVector.at(2 * i) - fullVector.at(2 * i + 1);
		v.push_back(s / 2);
	}

	return v;
}

std::vector<Mat_<float>> divideIntoFour(Mat_<float> originalImage)
{
	std::vector<Mat_<float>> results;
	int rows = originalImage.rows;
	int cols = originalImage.cols;

	Mat_<float> l = Mat_<float>(rows / 2, cols);
	Mat_<float> h = Mat_<float>(rows / 2, cols);

	for (int c = 0; c < cols; c++)
	{
		std::vector<float> half1Low;
		std::vector<float> half1High;
		for (int r = 0; r < rows; r++)
		{
			half1Low.push_back(originalImage(r, c));
			half1High.push_back(originalImage(r, c));
		}

		half1Low = getLowVector(half1Low);

		half1High = getHighVector(half1High);

		for (int r = 0; r < half1Low.size(); r++)
		{
			l(r, c) = half1Low.at(r);
			h(r, c) = half1High.at(r);
		}
	}

	Mat_<float> ll = Mat_<float>(rows / 2, cols / 2, 255);
	Mat_<float> lh = Mat_<float>(rows / 2, cols / 2, 255);
	Mat_<float> hl = Mat_<float>(rows / 2, cols / 2, 255);
	Mat_<float> hh = Mat_<float>(rows / 2, cols / 2, 255);

	for (int r = 0; r < rows / 2; r++)
	{
		std::vector<float> half1Low;
		std::vector<float> half1High;
		for (int c = 0; c < cols; c++)
		{
			half1Low.push_back(l(r, c));
			half1High.push_back(l(r, c));
		}

		half1Low = getLowVector(half1Low);

		half1High = getHighVector(half1High);


		for (int c = 0; c < half1Low.size(); c++)
		{
			ll(r, c) = half1Low.at(c);
			lh(r, c) = half1High.at(c);
		}
	}

	for (int r = 0; r < rows / 2; r++)
	{
		std::vector<float> half1Low;
		std::vector<float> half1High;
		for (int c = 0; c < cols; c++)
		{
			half1Low.push_back(h(r, c));
			half1High.push_back(h(r, c));
		}

		half1Low = getLowVector(half1Low);

		half1High = getHighVector(half1High);

		for (int c = 0; c < half1Low.size(); c++)
		{
			hl(r, c) = half1Low.at(c);
			hh(r, c) = half1High.at(c);
		}
	}

	results.push_back(ll);
	results.push_back(lh);
	results.push_back(hl);
	results.push_back(hh);

	return results;
}

std::vector<float> getLow_USample(std::vector<float> vector_low)
{
	std::vector<float> low_usample;
	int stop = 2 * vector_low.size();
	for (int k = 0; k < stop; k++)
	{
		low_usample.push_back(vector_low.at(k / 2));
	}
	return low_usample;
}

std::vector<float> getHigh_USample(std::vector<float> vector_high)
{
	std::vector<float> high_usample;
	int stop = 2 * vector_high.size();
	for (int k = 0; k < stop; k++)
	{
		high_usample.push_back(vector_high.at(k / 2) * hVec[k % 2]);
	}
	return high_usample;
}

Mat_<float> reconstructImage(Mat_<float> ll, Mat_<float> lh, Mat_<float> hl, Mat_<float> hh)
{
	int rows = ll.rows;
	int cols = ll.cols;
	Mat_<float> l = Mat_<float>(rows, cols * 2);
	Mat_<float> h = Mat_<float>(rows, cols * 2);
	for (int r = 0; r < rows; r++)
	{
		std::vector<float> vector_low;
		std::vector<float> vector_high;
		for (int c = 0; c < cols; c++)
		{
			vector_low.push_back(ll(r, c));
			vector_high.push_back(lh(r, c));
		}

		std::vector<float> low_usample = getLow_USample(vector_low);
		std::vector<float> high_usample = getHigh_USample(vector_high);
		for (int c = 0; c < low_usample.size(); c++)
		{
			l(r, c) = low_usample.at(c) + high_usample.at(c);
		}
	}

	for (int r = 0; r < rows; r++)
	{
		std::vector<float> vector_low;
		std::vector<float> vector_high;
		for (int c = 0; c < cols; c++)
		{
			vector_low.push_back(hl(r, c));
			vector_high.push_back(hh(r, c));
		}

		std::vector<float> low_usample = getLow_USample(vector_low);
		std::vector<float> high_usample = getHigh_USample(vector_high);

		for (int c = 0; c < low_usample.size(); c++)
		{
			h(r, c) = low_usample.at(c) + high_usample.at(c);
		}
	}

	Mat_<float> reconstructedImage = Mat_<float>(2 * rows, 2 * cols);
	for (int c = 0; c < 2 * cols; c++)
	{
		std::vector<float> vector_low;
		std::vector<float> vector_high;
		for (int r = 0; r < rows; r++)
		{
			vector_low.push_back(l(r, c));
			vector_high.push_back(h(r, c));
		}

		std::vector<float> low_usample = getLow_USample(vector_low);
		std::vector<float> high_usample = getHigh_USample(vector_high);

		for (int r = 0; r < low_usample.size(); r++)
		{
			reconstructedImage(r, c) = low_usample.at(r) + high_usample.at(r);
		}
	}

	return reconstructedImage;
}

// result = [LH_128x128, HL_128x128, HH_128x128, LH_64x64, HL_64X64, HH_64X64, ...., LL_2x2, LH_2X2, HL_2X2, HH_2X2]
// LL_2^n x 2^n -> LL_2^(n - 1), LH_2^(n - 1), HL_2^(n - 1), HH_2^(n - 1)
std::vector<Mat_<float>> recursiveDecomposition(Mat_<uchar> orig)
{
	std::vector<Mat_<float>> result;
	Mat_<uchar> ll = orig.clone();

	return result;
}

// allDecompositions = [LH_128x128, HL_128x128, HH_128x128, LH_64x64, HL_64X64, HH_64X64, ...., LL_2x2, LH_2X2, HL_2X2, HH_2X2]
// LL_4x4 = reconstructImage(LL_2x2, LH_2x2, HL_2x2, HH_2x2)
// ...
// LL_2^n = reconstructImage(LL_2^(n - 1), LH_2^(n - 1), HL_2^(n - 1), HH_2^(n - 1))
Mat_<float> recursiveReconstruction(std::vector<Mat_<float>> allDecompositions)
{
	Mat_<float> ll;
	return ll;
}

void recursiveTests()
{
	char fname[MAX_PATH];
	while (openFileDlg(fname))
	{
		Mat_<uchar> src = imread(fname, CV_LOAD_IMAGE_GRAYSCALE);
		std::vector<Mat_<float>> decomp = recursiveDecomposition(src);
		Mat_<uchar> reconstructed = recursiveReconstruction(decomp);
		imshow("Original", src);
		imshow("Reconstructed", reconstructed);
		waitKey(0);
	}
}

void testDecomposition()
{
	char fname[MAX_PATH];
	while (openFileDlg(fname))
	{
		Mat_<uchar> src = imread(fname, CV_LOAD_IMAGE_GRAYSCALE);
		std::vector<Mat_<float>> results = divideIntoFour(src);
		imshow("Original", src);

		Mat_<float> ll = results.at(0);
		Mat_<float> lh = results.at(1);
		Mat_<float> hl = results.at(2);
		Mat_<float> hh = results.at(3);
		Mat_<uchar> reconstructed = reconstructImage(ll, lh, hl, hh);

		Mat_<uchar> pll = ll;
		Mat_<uchar> plh = lh;
		Mat_<uchar> phl = hl;
		Mat_<uchar> phh = hh;

		imshow("Reconstruction", reconstructed);
		imshow("LL", pll);
		imshow("LH", plh + 128);
		imshow("HL", phl + 128);
		imshow("HH", phh + 128);

		waitKey(0);
	}
}

int main()
{
	testDecomposition();
	getchar();

	return 0;
}