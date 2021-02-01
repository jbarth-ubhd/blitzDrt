//#define BZ_DEBUG

#ifdef BZ_DEBUG
#warning BZ_DEBUG aktiviert!
#define DEBUG_MSG "BZ_DEBUG aktiviert! "
#else
#warning ohne BZ_DEBUG
#define DEBUG_MSG "ohne BZ_DEBUG: "
#endif

#include <cstdio>
#include <cmath>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <typeinfo>
#include <stdexcept>
#include <complex>
#include <algorithm>
#include <Magick++.h>
#include <fftw3.h>
#include <sys/time.h>
#include <boost/any.hpp>
#include <boost/program_options.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/numeric/conversion/bounds.hpp>
#include <boost/limits.hpp>
#include <blitz/array.h>
#include "fm3.h"
#include "sourcecode.h"

 template <class T>
 inline std::string stringify (const T& t)
 {
	  std::stringstream ss;
	  ss << t;
	  return ss.str();
 }

using namespace std;
using namespace Magick;
using namespace blitz;
using namespace boost::program_options;
using namespace boost::algorithm;

template<typename T> inline T sqr(T x) { return x*x; };
template<typename T> inline T identity(T x) { return x; }
inline double deg2rad(double x) { return x*M_PI/180.0; };
inline double rad2deg(double x) { return x*180.0/M_PI; };

double y_mxPlusB(double x, double x1, double y1, double x2, double y2) {
	double m=(y2-y1)/(x2-x1);
	double b=y1-m*x1;
	return m*x+b;
}

template<typename T1, typename T2> void checkBound(Array<T1, 2> &a, Array<T2, 2> &b) {
	for(int i=0; i<2; i++) {
		FIF(a.lbound(i)!=b.lbound(i));
		FIF(a.ubound(i)!=b.ubound(i));
	}
};

void imageCacheToArray(unsigned char *ic, Array<int, 2> &a) {
	for(int y=a.ubound(1); y>=a.lbound(1); y--) {
		for(int x=a.lbound(0); x<=a.ubound(0); x++) {
			a(x, y)=*(ic++);
		}
	}
}

void plotAsciiPixel(double v) {
	int x=int(v*10.999);
	if(x>10) { printf("^"); return; }
	if(x<0) { printf("v"); return; }
	switch(x) {
		case 10: printf(" "); break;
		case 9: printf("·"); break;
		case 8: printf(":"); break;
		case 7: printf("¬"); break;
		case 6: printf("÷"); break;
		case 5: printf("Y"); break;
		case 4: printf("I"); break;
		case 3: printf("U"); break;
		case 2: printf("D"); break;
		case 1: printf("N"); break;
		case 0: printf("Ø"); break;
	}
}

void printArrayAscii(Array<int, 2> &a) {
	printf("%d\n", a.ubound(1)-a.lbound(1));
	for(int y=a.ubound(1); y>=a.lbound(1)+1; y-=2) {
		printf("%5d:", y);
		for(int x=a.lbound(0); x<=a.ubound(0); x++) {
			plotAsciiPixel((a(x, y)+a(x, y-1))/2.0/255.0);
		}
		printf(" :%d\n", y);
	}
}

inline int intImageCenterMax(int x) { return (x-1)/2; }
inline int intImageCenterMin(int x) { return intImageCenterMax(x)-(x-1); }

void derivateTo(
		int (*diffTransform)(int),
		Array<int, 2> &from,
		Array<int, 2> &deri,
		double rxD, double ryD) {

	checkBound(from, deri);

	int rx=int(round(rxD));
	int ry=int(round(ryD));
	for(int x=from.lbound(0); x<=from.ubound(0); x++) {
		for(int y=from.lbound(1); y<=from.ubound(1); y++) {
			Array<int, 2> window(from,
				Range(max(from.lbound(0), x-rx),
					min(from.ubound(0), x+rx)),
				Range(max(from.lbound(1), y-ry),
					min(from.ubound(1), y+ry))
					); // Achtung: Konstruktor -> Referenz, NICHT Kopie!
			int m=int(round(mean(window)));
			deri(x, y)=diffTransform(from(x, y)-m);
		}
	}
}


class dRange: public Range {
public:
	int iMin, iMax;
	double dMin, dMax;

	dRange(int iMin, int iMax, double dMin, double dMax):
		Range(iMin, iMax),
		iMin(iMin), iMax(iMax),
		dMin(dMin), dMax(dMax) {}

	double indexToReal(double x) {
		return y_mxPlusB(x, iMin, dMin, iMax, dMax);
		//return (x-iMin)/(iMax-iMin)*(dMax-dMin)+dMin;
	}

	int realToIndex(double x) {
		return int(round(
				y_mxPlusB(x, dMin, iMin, dMax, iMax)
				));
		/*return int(round(
				(x-dMin)/(dMax-dMin)*(iMax-iMin)+iMin
				));*/
	}
};

int theoreticMaxDForImage(Array<int, 2> &a) {
	return int(ceil(sqrt(
			sqr(a.ubound(0))+sqr(a.ubound(1))
			)));
}

int theoreticMinDForImage(Array<int, 2> &a) {
	return int(floor(-sqrt(
			sqr(a.lbound(0))+sqr(a.lbound(1))
			)));
}

inline double sin90(double x) {	return sin(x+M_PI/2); }
inline double cos90(double x) { return cos(x+M_PI/2); }

void imageToDrt(Array<int, 2> &img, Array<int, 2> &drt, dRange &thetaRange) {
	for(int thetaIndex=drt.lbound(0); thetaIndex<=drt.ubound(0); thetaIndex++) {
		double theta=thetaRange.indexToReal(thetaIndex);

		Array<int, 1> summandYA(Range(img.lbound(1), img.ubound(1)));
		for(int y=img.lbound(1); y<=img.ubound(1); y++)
			summandYA(y)=int(round(65536.0*(sin90(theta)*y +0.5)));

		for(int x=img.lbound(0); x<=img.ubound(0); x++) {
			int summandX=int(round(65536.0*cos90(theta)*x));
			for(int y=img.lbound(1); y<=img.ubound(1); y++) {
				drt(thetaIndex, (summandYA(y)+summandX)>>16)
					+=img(x, y);
			}
		}

	}
}

void normalize(Array<int, 2> &a, int toMaxVal) {
	int aMin=min(a);
	int aMax=max(a);
	a=((a-aMin)*toMaxVal+(aMax-aMin)/2)/(aMax-aMin);
}

void normalize(Array<double, 2> &a) {
	double aMin=min(a);
	double aMax=max(a);
	a=(a-aMin)/(aMax-aMin);
}

void buildAxis(Array<int, 2> &drt, int t1, int d1, int t2, int d2, Array<int, 1> sharpDrtAxis) {
	double tStep=double(t2-t1)/(d2-d1);
	double t=t1+0.5;
	for(int d=d1; d<=d2; d++) {
		int t_int=int(floor(t));
		if(t_int<drt.lbound(firstDim)) t_int=drt.lbound(firstDim);
		if(t_int>drt.ubound(firstDim)) t_int=drt.ubound(firstDim);
		sharpDrtAxis(d)=drt(t_int, d);
		t+=tStep;
	}
}

void findMaxDrtAxis(Array<int, 2> &drt, dRange &thetaRange, double tStart, double tStop, Array<int, 2> &img, int *pT1, int *pT2, int *pT1GleichT2) {
	int indexTStart=thetaRange.realToIndex(tStart);
	int indexTStop =thetaRange.realToIndex(tStop);
	//!*!*!*!*!*!*! double Max=std::numeric_limits<double>::min(); // min()>0!!!
	double Max=boost::numeric::bounds<double>::lowest(); // -1E308
	double t1GleichT2Max=Max;
	for(int t1=indexTStart; t1<=indexTStop; t1++) {
		int d1=img.lbound(1);
		for(int t2=indexTStart; t2<=indexTStop; t2++) {
			int d2=img.ubound(1);
			Array<int, 1> sharpDrtAxis(Range(d1, d2));
			buildAxis(drt, t1, d1, t2, d2, sharpDrtAxis);
			double m=mean(sharpDrtAxis);
			if(m>=Max) { // *pT1+*pT2 werden sonst ggf. nie gesetzt - z. B. bei ganz weißem Bild.
				Max=m;
				*pT1=t1;
				*pT2=t2;
			}
			if(t1==t2 && m>=t1GleichT2Max) { // dito
				t1GleichT2Max=m;
				*pT1GleichT2=t1;
			}
		}
	}
}

void writeToImage(Array<double, 2> &a, Image &i) {
	for(int y=a.ubound(1), yImage=0; y>=a.lbound(1); y--, yImage++) {
		for(int x=a.lbound(0), xImage=0; x<=a.ubound(0); x++, xImage++) {
			if(a(x, y)<0) i.pixelColor(xImage, yImage, "blue");
			else
				if(a(x, y)>1) i.pixelColor(xImage, yImage, "red");
				else i.pixelColor(xImage, yImage, Magick::ColorGray(a(x, y)));
		}
	}
}

double calcSchwerpunktD(Array<int, 1> &sharpDrtAxis) {
	double sumD=0;
	double n=0;
	for(int d=sharpDrtAxis.lbound(0); d<=sharpDrtAxis.ubound(0); d++) {
		sumD+=d*sharpDrtAxis(d);
		n+=sharpDrtAxis(d);
	}
	return sumD/n;
}

class fftPoint {
public:
	double abs;
	double halfWavelength; /* nur historische Gründe */
	double lines;
	fftPoint() {}
	fftPoint(double a, double b, double c): abs(a), halfWavelength(b), lines(c) {}

};
bool operator<(const fftPoint& a, const fftPoint& b) {
    return a.abs < b.abs;
}

double sd(double sumX, double sumXX, double length) {
	if(length==0.0) { return 0.0; }
	return sqrt((sumXX-sqr(sumX)/length)/length);
}


class drtParams {
public:
	double zoomFactor;
	double graphW2;
	double graphH2;
	double graphW;
	double radT1;
	double radT2;
	int drtAxisFftLength;
	vector<fftPoint> drtAxisFftOrderedByAbsDesc;
	vector<double> drtAxis;
};

void getDrtParamsFromImage(variables_map &cmdLine, Image &graphImg, drtParams &drtParam) {
	double pmMaxTheta=deg2rad(cmdLine["angle"].as<double>());
	double thetaCompare=deg2rad(cmdLine["compareAngle"].as<double>());

	graphImg.type( GrayscaleType );

	graphImg.normalize();
	graphImg.negate(true);

	if (cmdLine.count("zoom")) {
		cout << "zooming image to " << cmdLine["zoom"].as<string>() << "\n";
		graphImg.scale(cmdLine["zoom"].as<string>());
	}
	drtParam.zoomFactor=graphImg.size().height()/(drtParam.graphH2*2);

	unsigned char *imageCache=new unsigned char[graphImg.size().width()*graphImg.size().height()];
	graphImg.write(0, 0, graphImg.size().width(), graphImg.size().height(), "R", CharPixel, imageCache);

	Range imageRangeX(intImageCenterMin(graphImg.size().width() ), intImageCenterMax(graphImg.size().width() ));
	Range imageRangeY(intImageCenterMin(graphImg.size().height()), intImageCenterMax(graphImg.size().height()));

	Array<int, 2> img(imageRangeX, imageRangeY);
	imageCacheToArray(imageCache, img);
	delete imageCache;

	Array<int, 2> sharpImg(imageRangeX, imageRangeY);

	printf("sharpening image..."); fflush(stdout);
	derivateTo(identity, img, sharpImg,
			graphImg.size().height()*cmdLine["sr"].as<double>()/100,
			graphImg.size().height()*cmdLine["sr"].as<double>()/100);
	printf("done.\n");

	double maxWinkelRad=pmMaxTheta +thetaCompare;
	dRange drtThetaRange(
			int(round(-maxWinkelRad*graphImg.size().height())),
			int(round( maxWinkelRad*graphImg.size().height())),
			-maxWinkelRad,
			 maxWinkelRad
			);

	dRange drtDRange(theoreticMinDForImage(sharpImg), theoreticMaxDForImage(sharpImg),
			-double(theoreticMinDForImage(sharpImg))/sharpImg.lbound(1)/2,
			double(theoreticMaxDForImage(sharpImg))/sharpImg.ubound(1)/2);

	Array<int, 2> drt(drtThetaRange, drtDRange);
	drt=0; // Initialisierung

	printf("imageToDrt..."); fflush(stdout);
	imageToDrt(sharpImg, drt, drtThetaRange);
	printf("done.\n");

	Array<int, 2> raiseDrt(drtThetaRange, drtDRange);
	printf("raising drt..."); fflush(stdout);
	derivateTo(abs, drt, raiseDrt,
			graphImg.size().height()*cmdLine["sr"].as<double>()/100,
			graphImg.size().height()*cmdLine["sr"].as<double>()/100);

	printf("sharpening drt..."); fflush(stdout);
	Array<int, 2> sharpDrt(drtThetaRange, drtDRange);
	derivateTo(identity, drt, sharpDrt,
			graphImg.size().height()*cmdLine["sr"].as<double>()/100,
			graphImg.size().height()*cmdLine["sr"].as<double>()/100);
	printf("done.\n");

	int t1, t2, t1GleichT2;
	int d1=sharpImg.lbound(1);
	int d2=sharpImg.ubound(1);
	printf("find max(raiseDrt)..."); fflush(stdout);

	findMaxDrtAxis(raiseDrt, drtThetaRange,
		-pmMaxTheta,
		 pmMaxTheta,
		sharpImg, &t1, &t2, &t1GleichT2);
	printf("done.\n");

	drtParam.radT1=drtThetaRange.indexToReal(t1);
	drtParam.radT2=drtThetaRange.indexToReal(t2);

	printf("angle at 0%%: %f; angle at 100%%: %f degrees\n", rad2deg(drtParam.radT1), rad2deg(drtParam.radT2));
	printf("max. var. at constant angle: %lf degrees\n", rad2deg(drtThetaRange.indexToReal(t1GleichT2)));

	//

	Array<int, 1> sharpDrtAxis(Range(d1, d2));
	buildAxis(sharpDrt, t1, d1, t2, d2, sharpDrtAxis);

	double doubleDrtAxis[d2-d1+1];
	printf("iDrtAxisVal:");
	for(int d=d1; d<=d2; d++) {
		doubleDrtAxis[d-d1]=double(sharpDrtAxis(d))/(sharpImg.ubound(0)-sharpImg.lbound(0)+1);
		printf(" %f", doubleDrtAxis[d-d1]);
		drtParam.drtAxis.push_back(doubleDrtAxis[d-d1]);
	}
	printf("\n");

	fftw_plan fftPlan;
	drtParam.drtAxisFftLength=(d2-d1+1)/2+1;
	complex<double>* fft = new complex<double>[drtParam.drtAxisFftLength];
	drtParam.drtAxisFftOrderedByAbsDesc=vector<fftPoint>(drtParam.drtAxisFftLength-1);
	//vector<fftPoint> drtParam.drtParam.drtAxisFftOrderedByAbsDesc(drtAxisFftLength-1);
/*  C++ has its own complex<T> template class, defined in the standard <complex> header file.
	  Reportedly, the C++ standards committee has recently agreed to mandate that the storage
	format used for this type be binary-compatible with the C99 type, i.e. an array T[2] with consecutive
	real [0] and imaginary [1] parts. (See report WG21/N1388.)
	  Although not part of the official standard as of this writing, the proposal stated that:
	“This solution has been tested with all current major implementations of the standard library and
	shown to be working.”
	  To the extent that this is true, if you have a variable complex<double> *x, you can pass it directly
	to FFTW via reinterpret_cast<fftw_complex*>(x).*/
	fftPlan = fftw_plan_dft_r2c_1d(d2-d1+1, doubleDrtAxis, reinterpret_cast<fftw_complex*>(fft), FFTW_ESTIMATE);
	fftw_execute(fftPlan);
	printf("iDrtAxisFFT(n=%d):", d2-d1+1);
	for(int i=1; i<drtParam.drtAxisFftLength; i++) {
		drtParam.drtAxisFftOrderedByAbsDesc[i-1].abs=abs(fft[i])/drtParam.drtAxisFftLength;
		// nächste Zeile: 33% (=0.33) war der default-Zoomfaktor der alten drt-Version.
		// damit die Zahlen vergleichbar sind, wird umgerechnet.
		drtParam.drtAxisFftOrderedByAbsDesc[i-1].halfWavelength=double(drtParam.drtAxisFftLength)*0.33/drtParam.zoomFactor/i;
		drtParam.drtAxisFftOrderedByAbsDesc[i-1].lines=i;
		printf(" %f", drtParam.drtAxisFftOrderedByAbsDesc[i-1].abs);
	}
	printf("\n");
	fftw_destroy_plan(fftPlan);
	delete[] fft;

	//

	Array<int, 1> raiseDrtAxis(Range(d1, d2));
	buildAxis(raiseDrt, t1, d1, t2, d2, raiseDrtAxis);
	double dSchwerRel=(calcSchwerpunktD(raiseDrtAxis)-d1)/(d2-d1);

	Array<int, 1> raiseDrtAxisVar1(Range(d1, d2));
	buildAxis(raiseDrt,
		t1-drtThetaRange.realToIndex(thetaCompare*   dSchwerRel ), d1,
		t2+drtThetaRange.realToIndex(thetaCompare*(1-dSchwerRel)), d2, raiseDrtAxisVar1);

	Array<int, 1> raiseDrtAxisVar2(Range(d1, d2));
	buildAxis(raiseDrt,
		t1+drtThetaRange.realToIndex(thetaCompare*   dSchwerRel ), d1,
		t2-drtThetaRange.realToIndex(thetaCompare*(1-dSchwerRel)), d2, raiseDrtAxisVar2);

	double var1qual=mean(raiseDrtAxisVar1)/mean(raiseDrtAxis);
	double var2qual=mean(raiseDrtAxisVar2)/mean(raiseDrtAxis);
	printf("drtAxis Var1:%f%% Var2:%f%% => %f\n",
			var1qual*100, var2qual*100, 100*(1-max(var1qual, var2qual)));

	if(cmdLine.count("drt")) {
		string sizeStr=stringify(drtThetaRange.iMax-drtThetaRange.iMin+1)
							+"x"+stringify(drtDRange.iMax-drtDRange.iMin+1);
		Image drtOut(sizeStr.c_str(), "black");
		drtOut.quality(100);

		Array<double, 2> tmpDrt(drtThetaRange, drtDRange);
		tmpDrt=raiseDrt*1.0;
		normalize(tmpDrt);

		double tStep=double(t2-t1)/(d2-d1);
		double t=t1+0.5;
		for(int d=d1; d<=d2; d++) {
			tmpDrt(int(floor(t)), d)=-1.0;
			t+=tStep;
		}

		writeToImage(tmpDrt, drtOut);

		try {
			drtOut.write(cmdLine["drt"].as<string>());
		}
		catch( Magick::WarningCoder &warning ) {
			cerr << "Coder Warning: " << warning.what() << endl;
		}
	    catch( Magick::Warning &warning ) {
	    	cerr << "Warning: " << warning.what() << endl;
	    }
	    catch( Magick::ErrorBlob &error) {
	    	cerr << "Error: " << error.what() << endl;
	    }
	}
}

string &after(string &haystack, const char *key, string &val) {
	if(haystack.find(key, 0)!=0) {
	    val=string("");
	    return val;
	}
	val=haystack.substr(strlen(key));
	return val;
}

int getDrtParamsFromDrtlog(variables_map &cmdLine, Image &graphImg, drtParams &drtParam) {
	drtParam.zoomFactor=0.33;

	ifstream drtLog;
	drtLog.open(cmdLine["import"].as<string>().c_str());
	if(drtLog.fail()) return 0;

	int paramcount=0;
	while(! drtLog.eof()) {
		string line;
		getline(drtLog, line);
		string val;

//			graphW/H wird nicht mehr angefasst, da bereits beim Laden des Bildes bestimmt.

		if((after(line, "size: ", val)).size()>0) {
			vector<string> size;
			split(size, val, is_any_of("x"));
			if(strtoul(size[0].c_str(), NULL, 10) != graphImg.size().width()
			|| strtoul(size[1].c_str(), NULL, 10) != graphImg.size().height()) {
				printf("FEHLER: Bildgröße aus %s stimmt nicht mit tatsächlicher Bildgröße überein.\n", cmdLine["import"].as<string>().c_str());
				return 0;
			}
		}


		if((after(line, "angle at ", val)).size()>0) {
			vector<string> angleParam;
			split(angleParam, val, is_any_of("%:; ")); // trennzeichen zusammenfassen

			double perD1=atof(angleParam[0].c_str());
			double radT1=deg2rad(atof(angleParam[3].c_str()));
			double perD2=atof(angleParam[7].c_str());
			double radT2=deg2rad(atof(angleParam[10].c_str()));

			drtParam.radT1=y_mxPlusB(  0, perD1, radT1, perD2, radT2);
			drtParam.radT2=y_mxPlusB(100, perD1, radT1, perD2, radT2);
			printf("radT1=%f\n", drtParam.radT1);
			printf("radT2=%f\n", drtParam.radT2);
			paramcount++;
		}

		if((after(line, "iDrtAxisFFT(n=", val)).size()>0) {
			vector<string> sfft;
			split(sfft, val, is_any_of(" "));
			vector<string> n;
			split(n, sfft[0], is_any_of(")"));
			int fftInputLength=atoi(n[0].c_str());
			drtParam.drtAxisFftLength=(fftInputLength-1)/2;
			// ja, war wohl in altem drt-Programm so;
			// die paar % machen den Kohl nicht fett.
			printf("drtParam.drtAxisFftLength=%d\n", drtParam.drtAxisFftLength);
			for(unsigned i=1; i<sfft.size(); i++) {
				drtParam.drtAxisFftOrderedByAbsDesc.push_back(fftPoint(
						atof(sfft[i].c_str()),
						double(drtParam.drtAxisFftLength)/i,
						double(i)
						));
			}
			paramcount++;
		}

		if((after(line, "iDrtAxisVal: ", val)).size()>0) {
			vector<string> string_drtAxis;
			split(string_drtAxis, val, is_any_of(" "));
			for(unsigned i=0; i<string_drtAxis.size(); i++) {
				drtParam.drtAxis.push_back(atof(string_drtAxis[i].c_str()));
			}
			paramcount++;
		}

		// "shrinking"="zooming"
		if((after(line, "shrinking image to ", val)).size()>0) {
			vector<string> zoom;
			split(zoom, val, is_any_of("%"));
			drtParam.zoomFactor=atof(zoom[0].c_str())/100;
			printf("zoomFactor=%f\n", drtParam.zoomFactor);
			paramcount++;
		}
		if((after(line, "zooming image to ", val)).size()>0) {
			vector<string> zoom;
			split(zoom, val, is_any_of("%"));
			drtParam.zoomFactor=atof(zoom[0].c_str())/100;
			printf("zoomFactor=%f\n", drtParam.zoomFactor);
			paramcount++;
		}
	}

	if(paramcount!=4) return 0;
	return 1;
}


int main(int argc, char **argv)
{
	FIF(sizeof(int)!=4);
	FIF((3>>1)!=1);
	FIF((-3>>1)!=-2);

	options_description optDesc(string(DEBUG_MSG)+string(__FILE__)+" @ "+string(__DATE__)+" "+string(__TIME__)+"\nAllowed options");
	optDesc.add_options()
	    ("help", "produce help message")
	    ("source", "print program source code")
	    ("if",  value<string>(), "input image file")
	    ("import", value<string>(), "import drtLog file")
	    ("drt", value<string>(), "drt image file (visualization of drt)")
	    ("tmpfile",  value<string>(), "temporary image file (format conversion only); will not be deleted")
	    ("of",  value<string>(), "output image file (after spherical correction)")
	    ("python", value<string>()->default_value("python"), "python command (requires PIL)")
	    ("doTrans", "does format conversion with transCmd, if applicable (sd<maxSd)")
	    ("noTransDeskew", "format conversion&statistics only")
	    ("maxSd", value<double>()->default_value(0.8), "max. sd(halfWavelength) allowing transformation")
	    ("sr",  value<double>()->default_value(1.0), "shapen radius [image height %]")
	    ("zoom",  value<string>()->default_value("33%"), "zoom +'%'")
	    ("angle",  value<double>()->default_value(5.0), "drt angle probe range [Grad]")
	    ("compareAngle",  value<double>()->default_value(1.0), "+/- angle variation to compare")
	    ("halfWavelengthVarianceAbsFactor",  value<double>()->default_value(50), "variance(halfWavelength(fft where abs(fft)>=max(abs(fft))*hWVAF))")

	    ("minPageHeight",  value<double>()->default_value(150.0), "page height minus border[mm]")
	    ("maxPageHeight",  value<double>()->default_value(600.0), "page height minus border[mm]")
	    ("minLineSpaceSize", value<double>()->default_value(10.0), "[pt] (1/72 inch = 0.35277mm)")
	    ("maxLineSpaceSize", value<double>()->default_value(16.0), "[pt] (1/72 inch = 0.35277mm)")
	    ("minLineSpaceRatio", value<double>()->default_value(0.7))
		("maxLineSpaceRatio", value<double>()->default_value(2.0))
		("minLineLike", value<double>()->default_value(26.664827), "%")
		("hysteresis", value<double>()->default_value(4.278268))
		("maxSdLineSpace", value<double>()->default_value(6.272757))
		("gamma", value<double>()->default_value(2.2), "sRGB gamma correction before(1/x)&after(x) transform")
	;

	variables_map cmdLine;
	store(parse_command_line(argc, argv, optDesc), cmdLine);
	notify(cmdLine);

	if (cmdLine.count("help")) {
	    cout << optDesc << "\n";
	    return 1;
	}
	if (cmdLine.count("source")) {
		cerr << SOURCECODE_TIMESTAMP << endl;
	    fwrite(SOURCECODE, 1, SOURCECODE_LENGTH, stdout);
	    return 1;
	}

	FIF(!cmdLine.count("if"));

	Image graphImg;
	printf("reading image..."); fflush(stdout);
	try {
		graphImg.read(cmdLine["if"].as<string>()); printf("done.\n");
	}
	catch( Magick::WarningCoder &warning ) {
		cerr << "Coder Warning: " << warning.what() << endl;
	}
	catch( Magick::Warning &warning ) {
		cerr << "Warning: " << warning.what() << endl;
	}
	catch( Magick::ErrorBlob &error) {
		cerr << "Error: " << error.what() << endl;
	}

	Image transImg=graphImg; // wird unten gebraucht

	graphImg.modifyImage();
	//SetImageVirtualPixelMethod( graphImg.image(), MagickCore::WhiteVirtualPixelMethod);
	//TODO: will in aktueller Version noch nicht
	// graphImg.options()->virtualPixelMethod( MagickCore::WhiteVirtualPixelMethod );

	graphImg.virtualPixelMethod(MagickCore::WhiteVirtualPixelMethod); // aber auch das will in Magick++ unter (ubuntu) karmic koala nicht

	drtParams drtParam;
	printf("size: %lux%lu\n", graphImg.size().width(), graphImg.size().height());
	drtParam.graphW2=double(graphImg.size().width())/2;
	drtParam.graphH2=double(graphImg.size().height())/2;
	drtParam.graphW=drtParam.graphW2*2;

	if(cmdLine.count("tmpfile")) {
		// ist ja nur ein tmpfile: graphImg.density(Geometry(300, 300)); // 300 dpi setzen
		graphImg.quality(100);
		printf("saving image..."); fflush(stdout);
		//graphImg.write(cmdLine["tmpfile"].as<string>());
		printf("done.\n");
	}

	if(!cmdLine.count("import")) {
		getDrtParamsFromImage(cmdLine, graphImg, drtParam);
	} else {
		if(!getDrtParamsFromDrtlog(cmdLine, graphImg, drtParam)) {
			fprintf(stderr, "WARNUNG: drtLog-(import)-Datei %s fehlerhaft. DRT-Analyse wird vom Bild durchgefuehrt.\n", cmdLine["import"].as<string>().c_str());
			getDrtParamsFromImage(cmdLine, graphImg, drtParam);
		}
	}

	double radT12avg=(drtParam.radT1+drtParam.radT2)/2;
	double radT1rel=drtParam.radT1-radT12avg;
	double radT2rel=drtParam.radT2-radT12avg;

	// 0 3
	// 1 2

	vector<complex<double> > transKoord(4);

	if(radT2rel<0) {
		transKoord[0]=complex<double>(-drtParam.graphW2, drtParam.graphH2);
		transKoord[3]=transKoord[0]+polar(drtParam.graphW, radT2rel);
	} else {
		transKoord[3]=complex<double>(drtParam.graphW2, drtParam.graphH2);
		transKoord[0]=transKoord[3]-polar(drtParam.graphW, radT2rel);
	}

	if(radT1rel>0) {
		transKoord[1]=complex<double>(-drtParam.graphW2, -drtParam.graphH2);
		transKoord[2]=transKoord[1]+polar(drtParam.graphW, radT1rel);
	} else {
		transKoord[2]=complex<double>(drtParam.graphW2, -drtParam.graphH2);
		transKoord[1]=transKoord[2]-polar(drtParam.graphW, radT1rel);
	}

	double distortParam[16];
	for(int i=0; i<4; i++) {
		int x, y;
		switch(i) {
			case 0: x=0; y=0; break;
			case 1: x=0; y=drtParam.graphH2*2-1; break;
			case 2: x=drtParam.graphW2*2-1; y=drtParam.graphH2*2-1; break;
			case 3: x=drtParam.graphW2*2-1; y=0; break;
		}
		double r=abs(transKoord[i]);
		double phi=std::arg(transKoord[i])+radT12avg;
		distortParam[i*4+0]=r*cos(phi)+drtParam.graphW2; //+pilXoff;
		distortParam[i*4+1]=(drtParam.graphH2*2-1)-(r*sin(phi)+drtParam.graphH2); //+pilYoff;
		distortParam[i*4+2]=x;
		distortParam[i*4+3]=y;

		printf("%d: %f,%f %f,%f\n", i, distortParam[i*4+0], distortParam[i*4+1], distortParam[i*4+2], distortParam[i*4+3]);
	}

	sort(drtParam.drtAxisFftOrderedByAbsDesc.begin(), drtParam.drtAxisFftOrderedByAbsDesc.end());
	reverse(drtParam.drtAxisFftOrderedByAbsDesc.begin(), drtParam.drtAxisFftOrderedByAbsDesc.end());
	printf("iDrtAxisFFT ordered by abs: ");
	double sumW=0;
	double sumL=0;
	double sumWW=0;
	int n=0;
	for(int i=0; i<drtParam.drtAxisFftLength-1
		&& drtParam.drtAxisFftOrderedByAbsDesc[i].abs>=
			drtParam.drtAxisFftOrderedByAbsDesc[0].abs
				*cmdLine["halfWavelengthVarianceAbsFactor"].as<double>()/100;
		i++) {
		sumW+=drtParam.drtAxisFftOrderedByAbsDesc[i].halfWavelength;
		sumWW+=pow(drtParam.drtAxisFftOrderedByAbsDesc[i].halfWavelength, 2);
		sumL+=drtParam.drtAxisFftOrderedByAbsDesc[i].lines;
		n++;
		printf(" (%.3lf:%.3lf)", drtParam.drtAxisFftOrderedByAbsDesc[i].abs,
				drtParam.drtAxisFftOrderedByAbsDesc[i].halfWavelength);
	}
	printf("\n");
	printf("sd=%lf halfWl=%lf lines=%lf\n", sd(sumW, sumWW, n), sumW/n, sumL/n);

	// neue Evaluierung der Drt-Achse:

	double minLines=cmdLine["minPageHeight"].as<double>() / (cmdLine["maxLineSpaceSize"].as<double>()*0.35277);
	double maxLines=cmdLine["maxPageHeight"].as<double>() / (cmdLine["minLineSpaceSize"].as<double>()*0.35277);
	int minLineSpace=int(drtParam.drtAxis.size()/maxLines);
	int maxLineSpace=int(drtParam.drtAxis.size()/minLines+0.999);
	printf("minLineSpace:%d\n", minLineSpace);
	printf("maxLineSpace:%d\n", maxLineSpace);
	double minRatio=cmdLine["minLineSpaceRatio"].as<double>();
	double maxRatio=cmdLine["maxLineSpaceRatio"].as<double>();

	double line=0;
	double space=0;
	double sumLine=0;
	double sumLine2=0;
	double sumSpace=0;
	double sumSpace2=0;
	double sumN=0;

	int state=0;
	double hyst=cmdLine["hysteresis"].as<double>();

	double insgesamt=0;
	double zeilenartig=0;

	if(drtParam.drtAxis[0]>0) { state=1; space++; } else { state=0; line++; }

	printf("iDrtAxis-LineSpacePattern:");
	for(unsigned i=1; i<drtParam.drtAxis.size(); i++) {
		if(state==0) {
			if(drtParam.drtAxis[i]>hyst) {
				printf("(l:%.0lf;s:%.0lf)s", line, space);
				if(space==0.0) { space+=0.000001; }
				double ratio=line/space;
				insgesamt+=line+space;
				if(line+space>=minLineSpace && line+space<=maxLineSpace)	{
					if(ratio>=minRatio/1.55 && ratio<=maxRatio*1.55) {
						sumLine+=line;
						sumLine2+=line*line;
						sumSpace+=space;
						sumSpace2+=space*space;
						sumN++;
					}
					if(ratio>=minRatio && ratio<=maxRatio) {
						zeilenartig+=line+space;
					}
				}
				line=0; space=1;
				state=1;
			} else { printf("l"); line++; } // hyst. n. überschr.
		} else { // status=1
			if(drtParam.drtAxis[i]<(-hyst)) {
				line++;
				printf("l");
				state=0;
			} else { printf("s"); space++; }
		}
	}
	printf("\n");
	double sdLine=sd(sumLine, sumLine2, sumN);
	double sdSpace=sd(sumSpace, sumSpace2, sumN);
	double sdLineSpace=sqrt(sdLine*sdLine+sdSpace*sdSpace);

	double proz_zeilenartig=double(zeilenartig)*100.0/(insgesamt+boost::numeric::bounds<double>::smallest()); // sicherheitshalber, falls weißes Blatt etc.
	printf("line-like: %.2f%%; sdLineSpace=%.3f\n", proz_zeilenartig, sdLineSpace);

	if(cmdLine.count("doTrans") && cmdLine.count("of")) {
		//if(sd(sumW, sumWW, n) < cmdLine["maxSd"].as<double>() && cmdLine.count("noTransDeskew")==0)
		// wird nun in doUnzip...sh erledigt: transImg.density(Geometry(300, 300)); // 300 dpi setzen
		if(proz_zeilenartig >= cmdLine["minLineLike"].as<double>()
				&& sdLineSpace <= cmdLine["maxSdLineSpace"].as<double>()
				&& cmdLine.count("noTransDeskew")==0) {
			transImg.compressType(NoCompression);
			transImg.depth(16);
			transImg.gamma(1/cmdLine["gamma"].as<double>());
			transImg.distort(MagickCore::PerspectiveDistortion,  16, distortParam);
			transImg.gamma(cmdLine["gamma"].as<double>());
			transImg.depth(8);
		}
		try {
			transImg.write(cmdLine["of"].as<string>());
		}
		catch( Magick::WarningCoder &warning ) {
			cerr << "Coder Warning: " << warning.what() << endl;
		}
		catch( Magick::Warning &warning ) {
			cerr << "Warning: " << warning.what() << endl;
		}
		catch( Magick::ErrorBlob &error) {
			cerr << "Error: " << error.what() << endl;
		}
	}
}

