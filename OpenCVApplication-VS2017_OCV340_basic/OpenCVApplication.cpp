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

										   // the fastest approach using the “diblook style”
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

		std::vector<float> low_usample;
		std::vector<float> high_usample;

		int stop = 2 * vector_low.size();
		for (int k = 0; k < stop; k++)
		{
			low_usample.push_back(vector_low.at(k / 2));
			high_usample.push_back(vector_high.at(k / 2) * hVec[k % 2]);
		}
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

		std::vector<float> low_usample;
		std::vector<float> high_usample;

		int stop = 2 * vector_low.size();
		for (int k = 0; k < stop; k++)
		{
			low_usample.push_back(vector_low.at(k / 2));
			high_usample.push_back(vector_high.at(k / 2) * hVec[k % 2]);
		}
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

		std::vector<float> low_usample;
		std::vector<float> high_usample;

		int stop = 2 * vector_low.size();
		for (int k = 0; k < stop; k++)
		{
			low_usample.push_back(vector_low.at(k / 2));
			high_usample.push_back(vector_high.at(k / 2) * hVec[k % 2]);
		}
		for (int r = 0; r < low_usample.size(); r++)
		{
			reconstructedImage(r, c) = low_usample.at(r) + high_usample.at(r);
		}
	}

	return reconstructedImage;
}

std::vector<Mat_<float>> recursiveDecomposition(Mat_<uchar> orig)
{
	std::vector<Mat_<float>> result;
	Mat_<uchar> ll = orig.clone();
	while (ll.rows > 2)
	{
		std::vector<Mat_<float>> div4 = divideIntoFour(ll);
		result.push_back(div4.at(0));
		
		result.push_back(div4.at(1));
		result.push_back(div4.at(2));
		result.push_back(div4.at(3));
		ll = div4.at(0).clone();
	}

	return result;
}

// LL, LH, HL, HH apar cate 4, una dupa alta, pentru un anumit nivel de adancime
Mat_<float> recursiveReconstruction(std::vector<Mat_<float>> allDecompositions)
{
	int groupsOf4Count = allDecompositions.size() / 4;
	Mat_<float> ll = allDecompositions.at(allDecompositions.size() - 4).clone();
	for (int i = 0; i < groupsOf4Count; i++)
	{
		int start = allDecompositions.size() - 4 * (i + 1);
		Mat_<float> lh = allDecompositions.at(start + 1);
		Mat_<float> hl = allDecompositions.at(start + 2);
		Mat_<float> hh = allDecompositions.at(start + 3);
		Mat_<float> newLL = reconstructImage(ll, lh, hl, hh);
		ll = newLL;
	}
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
	//recursiveTests();
	getchar();

	return 0;
}