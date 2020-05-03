// OpenCVApplication.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "common.h"
#include <queue>
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



Mat_<uchar> combineImage(Mat_<uchar> ll, Mat_<uchar> lh, Mat_<uchar> hl, Mat_<uchar> hh) {

	Mat_<uchar> result(ll.rows * 2, ll.cols * 2);
	ll.copyTo(result(Rect(0, 0, ll.cols, ll.rows)));
	lh.copyTo(result(Rect(lh.rows, 0, lh.cols, lh.rows)));
	hl.copyTo(result(Rect(0, hl.rows, hl.cols, hl.rows)));
	hh.copyTo(result(Rect(hh.rows, hh.cols, hh.cols, hh.rows)));

	return result;
}
// result = [LH_128x128, HL_128x128, HH_128x128, LH_64x64, HL_64X64, HH_64X64, ...., LH_2X2, HL_2X2, HH_2X2, LL_2x2]
// LL_2^n x 2^n -> LL_2^(n - 1), LH_2^(n - 1), HL_2^(n - 1), HH_2^(n - 1)
std::vector<Mat_<float>> recursiveDecomposition(Mat_<uchar> orig)
{
	std::vector<Mat_<float>> result;
	Mat_<uchar> ll = orig.clone();

	while (ll.rows > 2) {
		std::vector<Mat_<float>> divFour = divideIntoFour(ll);
		ll = divFour.at(0).clone();
		result.push_back(divFour.at(1));
		result.push_back(divFour.at(2));
		result.push_back(divFour.at(3));
		if (ll.rows == 2) {
			result.push_back(ll);
		}
	}

	return result;
}

// allDecompositions = [LH_128x128, HL_128x128, HH_128x128, LH_64x64, HL_64X64, HH_64X64, ...., LH_2X2, HL_2X2, HH_2X2, LL2x2] -> LL2x2 sa fie ultimul
// LL_4x4 = reconstructImage(LL_2x2, LH_2x2, HL_2x2, HH_2x2)
// ...
// LL_2^n = reconstructImage(LL_2^(n - 1), LH_2^(n - 1), HL_2^(n - 1), HH_2^(n - 1))
Mat_<float> recursiveReconstruction(std::vector<Mat_<float>> allDecompositions)
{

	int sz = allDecompositions.size();
	Mat_<float> ll = reconstructImage(allDecompositions.at(sz - 1), allDecompositions.at(sz - 4), allDecompositions.at(sz - 3), allDecompositions.at(sz - 2));

	for (int i = sz - 5; i >= 0; i -= 3) {
		ll = reconstructImage(ll, allDecompositions.at(i - 2), allDecompositions.at(i - 1), allDecompositions.at(i));
	}

	return ll;
}


Mat_<uchar> modifyContrast(Mat_<uchar> img) {
	int imin = 0xffffff, imax = -0xffffff, outMax = 250, outMin = 10;

	for (int i = 0; i < img.rows; i++) {
		for (int j = 0; j < img.cols; j++) {
			//g in max/g in min
			if (img(i, j) < imin) {
				imin = img(i, j);
			}
			if (img(i, j) > imax) {
				imax = img(i, j);
			}
		}
	}
	float decision = (float)(outMax - outMin) / (float)(imax - imin);
	Mat_<uchar>contrast(img.rows, img.cols);
	for (int i = 0; i < img.rows; i++) {
		for (int j = 0; j < img.cols; j++) {
			contrast(i, j) = outMin + (img(i, j) - imin)*decision;
		}
	}
	return contrast;

}

// result = [LL_128x128,LH_128x128, HL_128x128, HH_128x128,LL_64x64, LH_64x64, HL_64X64, HH_64X64, ...., LL_16x16, LH_16X16, HL_16X16, HH_16X16]

void display4Levels() {
	char fname[MAX_PATH];
	while (openFileDlg(fname))
	{
		Mat_<uchar> src = imread(fname, CV_LOAD_IMAGE_GRAYSCALE);
		std::vector<Mat_<uchar>> finalImg;
		imshow("Original", src);
		for (int i = 0; i < 4; i++) {
			std::vector<Mat_<float>> results = divideIntoFour(src);

			Mat_<float> ll = results.at(0);
			Mat_<float> lh = results.at(1);
			Mat_<float> hl = results.at(2);
			Mat_<float> hh = results.at(3);


			Mat_<uchar> pll = ll;
			Mat_<uchar> plh = lh;
			Mat_<uchar> phl = hl;
			Mat_<uchar> phh = hh;

			finalImg.push_back(pll);
			finalImg.push_back(plh);
			finalImg.push_back(phl);
			finalImg.push_back(phh);

			src = pll;
		}
		Mat_<uchar> combineImg;
		for (int i = finalImg.size() - 1; i >= 0; i -= 4) {
			Mat_<uchar> pll = finalImg.at(i - 3);
			Mat_<uchar> plh = finalImg.at(i - 2);
			Mat_<uchar> phl = finalImg.at(i - 1);
			Mat_<uchar> phh = finalImg.at(i);
			if (pll.rows == 16) {
				combineImg = combineImage(pll, modifyContrast(plh), modifyContrast(phl), modifyContrast(phh));
			}
			else {
				combineImg = combineImage(combineImg, modifyContrast(plh), modifyContrast(phl), modifyContrast(phh));
			}

		}
		imshow("256x256", combineImg);
		waitKey(0);
	}
}

void testRecursiveReconstruction() {
	char fname[MAX_PATH];
	while (openFileDlg(fname))
	{
		Mat_<uchar> src = imread(fname, CV_LOAD_IMAGE_GRAYSCALE);
		std::vector<Mat_<float>> decompositions = recursiveDecomposition(src);
		Mat_<float> finalImg = recursiveReconstruction(decompositions);
		Mat_<uchar> reconstructed = finalImg;

		imshow("Original", src);
		imshow("Reconstruction", reconstructed);
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
		imshow("LH", modifyContrast(plh));
		imshow("HL", modifyContrast(phl));
		imshow("HH", modifyContrast(phh));

		waitKey(0);
	}
}

// res -> (original - reconstruction) * 10 + 128
Mat_<uchar> computeDifference(Mat_<uchar> original, Mat_<uchar> reconstruction)
{
	int rows = original.rows;
	int cols = original.cols;
	Mat_<uchar> res = Mat_<uchar>(rows, cols);
	for (int i = 0; i < original.rows; i++) {
		for (int j = 0; j < original.cols; j++) {
			res(i, j) = (original(i, j) - reconstruction(i, j)) * 10 + 128;
		}
	}
	//res = res * 10 + 128;
	return res;
}

// se afiseaza src (original), imaginea reconstruita (dupa divizare recursiva si reconstructie recursiva)
//		si imaginea diferenta returnata de functia de mai sus
void testOriginalComparisonWithRes()
{
	char fname[MAX_PATH];
	while (openFileDlg(fname))
	{
		Mat_<uchar> img = imread(fname, CV_LOAD_IMAGE_GRAYSCALE);
		std::vector<Mat_<float>> decompositions = recursiveDecomposition(img);
		Mat_<float> finalImg = recursiveReconstruction(decompositions);
		Mat_<uchar> reconstructed = finalImg;
		Mat_<uchar> dif = computeDifference(img, reconstructed);

		imshow("Original", img);
		imshow("Reconstruction", reconstructed);
		imshow("Difference", dif);

		waitKey(0);
	}
}

void printVector(std::vector<float> vec)
{
	std::cout << "[";
	for (int i = 0; i < vec.size() - 1; i++)
	{
		std::cout << vec.at(i) << ", ";
	}
	std::cout << vec.at(vec.size() - 1) << "]";
}

std::vector<float> addVectors(std::vector<float> v1, std::vector<float> v2)
{
	std::vector<float> res;
	for (int i = 0; i < v1.size(); i++)
	{
		res.push_back(v1.at(i) + v2.at(i));
	}
	return res;
}

// param nou -> nr nivele descomp
// nivel 0 -> nivel n-2 adaugare in vec de matrici rezultat doar high
// nivel n - 1 => adaugare in vec de matrici rezultat high & low
std::vector<std::vector<float>> printAndReturnDecompositionOfVec(std::vector<float> nums, int levels)
{
	std::vector<std::vector<float>> result;
	std::vector<float> lowVec = nums;
	std::vector<float> highVec;
	for (int i = 1; i <= levels; i++)
	{
		printVector(lowVec);
		std::cout << " -> ";
		highVec = getHighVector(lowVec);
		lowVec = getLowVector(lowVec);
		std::cout << "L : ";
		printVector(lowVec);
		std::cout << ", H: ";
		printVector(highVec);

		result.push_back(highVec);
		if (i == levels)
		{
			result.push_back(lowVec);
		}
		std::cout << "\n";
	}

	return result;
}

void reconstructTestVector(std::vector<std::vector<float>> decomp)
{
	int size = decomp.size();
	std::vector<float> low = decomp.at(size - 1);
	std::vector<float> high = decomp.at(size - 2);

	std::vector<float> lowUSample = getLow_USample(low);
	std::vector<float> highUSample = getHigh_USample(high);
	low = addVectors(lowUSample, highUSample);

	std::cout << "LowUS: ";
	printVector(lowUSample);
	std::cout << " + HighUS: ";
	printVector(highUSample);
	std::cout << " = Low: ";
	printVector(low);
	std::cout << std::endl;

	for (int i = size - 3; i >= 0; i--)
	{
		high = decomp.at(i);
		lowUSample = getLow_USample(low);
		highUSample = getHigh_USample(high);
		low = addVectors(lowUSample, highUSample);
		std::cout << "LowUS: ";
		printVector(lowUSample);
		std::cout << " + HighUS: ";
		printVector(highUSample);
		std::cout << " = Low: ";
		printVector(low);
		std::cout << std::endl;
	}

	std::cout << std::endl;
}

void testSimpleVector()
{
	std::vector<float> nums = { 9, 7, 3, 5, 6, 10, 2, 6 };
	std::vector<std::vector<float>> decomposition = printAndReturnDecompositionOfVec(nums, 3);
	std::cout << "\n";
	reconstructTestVector(decomposition);
}


// returneaza o noua matrice dst cu valorile:
	// dst(i, j) = 0 daca abs(src(i, j)) < threshold; dst(i, j) = src(i, j) altfel
Mat_<float> filterMatrixWithThreshold(Mat_<float> src, int threshold)
{
	return Mat_<float>();
}


//	decomposition = [LH_128x128, HL_128x128, HH_128x128, LH_64x64, HL_64X64, HH_64X64, ...., LH_2X2, HL_2X2, HH_2X2, LL2x2]
//	se returneaza un nou vector cu matrici noi (nu se modifica nici vectorul decomposition, nici vreo matrice din el)
//		in care matricile LH, HL si HH din decomposition(la orice nivel) se transforma cu functia de mai sus "filterMatrixWithThreshold"
//			(deci se adauga in vectorul rezultat o noua matrice transformata),
//		iar cand se ajunge la matricea LL_2x2, se adauga o copie a acesteia in vectorul rezultat

std::vector<Mat_<float>> filterHMatricesWithThreshold(std::vector<Mat_<float>> decomposition, int threshold)
{
	return std::vector<Mat_<float>>();
}


// Se incarca imagini in mod obisnuit
//		Imaginea incarcata se va descompune cu functia recursiveDecomposition
//		Se apeleaza functia filterHMatricesWithThreshold cu rezultatul acestei descompuneri si un prag citit de la tastatura
//			Rezultatul apelului functiei va fi transmis lui recursiveReconstruction
//				Se va afisa imaginea originala si imaginea rezultata in urma apelului recursiveReconstruction
void testNoiseFilter()
{

}

// Se completeaza meniul cu fiecare noua functionalitate
int main()
{
	int op;

	do {
		destroyAllWindows();
		printf("Menu:\n");
		printf(" 1 - Decomposition & Reconstruction \n");
		printf(" 2 - Recursive Reconstruction\n");
		printf(" 3 - Recursive 4 Levels Decomposition \n");
		printf(" 4 - Compare original to reconstructed\n");
		printf(" 5 - Test vector\n");
		printf(" 0 - Exit\n\n");
		printf("Option: ");
		scanf("%d", &op);

		switch (op) {
			case 1:
				testDecomposition();
				break;
			case 2:
				testRecursiveReconstruction();
				break;
			case 3:
				display4Levels();
				break;
			case 4:
				testOriginalComparisonWithRes();
				break;

			case 5:
				testSimpleVector();
				break;
		}
	} while (op != 0);
	
	getchar();

	return 0;
}