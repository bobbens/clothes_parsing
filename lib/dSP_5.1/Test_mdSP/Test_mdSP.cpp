/*  This file is part of distributed Structured Prediction (dSP) - http://www.alexander-schwing.de/
 *
 *  distributed Structured Prediction (dSP) is free software: you can
 *  redistribute it and/or modify it under the terms of the GNU General
 *  Public License as published by the Free Software Foundation, either
 *  version 3 of the License, or (at your option) any later version.
 *
 *  distributed Structured Prediction (dSP) is distributed in the hope
 *  that it will be useful, but WITHOUT ANY WARRANTY; without even the
 *  implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 *  PURPOSE. See the GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with distributed Structured Prediction (dSP).
 *  If not, see <http://www.gnu.org/licenses/>.
 *
 *  Copyright (C) 2010-2013  Alexander G. Schwing  [http://www.alexander-schwing.de/]
 */

//Author: Alexander G. Schwing

/*  ON A PERSONAL NOTE: I spent a significant amount of time to go through
 *  both, theoretical justifications and coding of this framework.
 *  I hope the package is useful for your task. Any requests, 
 *  feedback, donations and support would be greatly appreciated.
 *  Thank you for contacting me!
 */
#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <string.h>

#include "../libRegionBPP/libRegionBPP.h"
#include "../libDHMultM/libDHSampPMultM.h"
#include "../libDHMultM/libDHModlPMultM.h"
#include "../libSPRegion/libSPRegionMultM.h"

//#define WITH_FLOAT			//if defined: using float rather than double
//#define WITH_PARALLELINF		//if defined: parallel inference rather than parallel in the samples


#ifdef WITH_FLOAT
typedef float NUM_TYPE;
#else
typedef double NUM_TYPE;
#endif

#ifdef WITH_PARALLELINF
typedef libRegionBPP<NUM_TYPE,true> SOLVER_TYPE;
typedef libDHModlPMultM<NUM_TYPE,SOLVER_TYPE> DH_TYPE;
#else
typedef libRegionBPP<NUM_TYPE,false> SOLVER_TYPE;
typedef libDHSampPMultM<NUM_TYPE,SOLVER_TYPE> DH_TYPE;
#endif

typedef libSPRegionMultM<NUM_TYPE,libDHBaseMultM<NUM_TYPE> > SP_TYPE;

#ifdef USE_ON_WINDOWS
#include <io.h>
void TraverseDirectory(const std::string& path, std::string& pattern, bool subdirectories, std::vector<std::string>& fileNames) {
	struct _finddatai64_t data;
	std::string fname = path + "\\" + pattern;
	// start the finder -- on error _findfirsti64() will return -1, otherwise if no
	// error it returns a handle greater than -1.
	intptr_t h = _findfirsti64(fname.c_str(),&data);
	if(h >= 0) {
		do {
			if( (data.attrib & _A_SUBDIR) ) {
				if( subdirectories && strcmp(data.name,".") != 0 && strcmp(data.name,"..") != 0) {
					fname = path + "\\" + data.name;
					TraverseDirectory(fname,pattern,true, fileNames);
				}
			} else {
				fileNames.push_back(path + "\\" + data.name);
			}
		} while( _findnexti64(h,&data) == 0);

		_findclose(h);
	}
}
#else
#include <dirent.h>
#include <fnmatch.h>
void TraverseDirectory(const std::string& path, std::string& pattern, bool subdirectories, std::vector<std::string>& fileNames) {
	DIR *dir, *tstdp;
    struct dirent *dp;

    //open the directory
    if((dir  = opendir(path.c_str())) == NULL)
    {
		std::cout << "Error opening " << path << std::endl;
        return;
    }

    while ((dp = readdir(dir)) != NULL)
    {
        tstdp=opendir(dp->d_name);
		
		if(tstdp) {
			closedir(tstdp);
			if(subdirectories) {
				//TraverseDirectory(
			}
		} else {
			if(fnmatch(pattern.c_str(), dp->d_name, 0)==0) {
				std::string tmp(path);
				tmp.append("/").append(dp->d_name);
				fileNames.push_back(tmp);
				//std::cout << fileNames.back() << std::endl;
			}
		}
    }

    closedir(dir);
    return;

}
#endif

int ParseInput(int argc, char** argv, struct DH_TYPE::SPParameters& params, std::string& GTPath, int& task, char*& weightFile, char*& leaveFile) {
	for(int k=1;k<argc;++k) {
		if(strcmp(argv[k], "-c")==0 && k+1!=argc) {
			params.C = NUM_TYPE(atof(argv[++k]));
		} else if(::strcmp(argv[k], "-e")==0 && k+1!=argc) {
			params.epsilon = NUM_TYPE(atof(argv[++k]));
		} else if(::strcmp(argv[k], "-p")==0 && k+1!=argc) {
			params.p = NUM_TYPE(atof(argv[++k]));
		} else if(::strcmp(argv[k], "-a")==0 && k+1!=argc) {
			params.ArmijoIterations = size_t(atoi(argv[++k]));
		} else if(::strcmp(argv[k], "-j")==0 && k+1!=argc) {
			params.CRFIterations = size_t(atoi(argv[++k]));
		} else if(::strcmp(argv[k], "-k")==0 && k+1!=argc) {
			params.CRFMPIterations = size_t(atoi(argv[++k]));
		} else if(::strcmp(argv[k], "-f")==0 && k+1!=argc) {
			params.CRFEraseMessages = size_t(atoi(argv[++k]));
		} else if(::strcmp(argv[k], "-r")==0 && k+1!=argc) {
			params.ReuseMessagesForF = size_t(atoi(argv[++k]));
		} else if(::strcmp(argv[k], "-o")==0 && k+1!=argc) {
			params.CRFOuterIterations = size_t(atoi(argv[++k]));
		} else if(::strcmp(argv[k], "-x")==0 && k+1!=argc) {
			params.CRFOuterExchange = size_t(atoi(argv[++k]));
		} else if(::strcmp(argv[k], "-m")==0 && k+1!=argc) {
			params.MPIterations = size_t(atoi(argv[++k]));
		} else if(::strcmp(argv[k], "-b")==0 && k+1!=argc) {
			params.BetheCountingNumbers = atoi(argv[++k]);
		} else if(::strcmp(argv[k], "-i")==0 && k+1!=argc) {
			params.CCCPIterations = size_t(atoi(argv[++k]));
		} else if(::strcmp(argv[k], "-d")==0 && k+1!=argc) {
			GTPath = std::string(argv[++k]);
		} else if(::strcmp(argv[k], "-t")==0 && k+1!=argc) {
			task = atoi(argv[++k]);
		} else if(::strcmp(argv[k], "-y")==0 && k+1!=argc) {
			params.ReadBinary = atoi(argv[++k]);
		} else if(::strcmp(argv[k], "-w")==0 && k+1!=argc) {
			weightFile = argv[++k];
		} else if(::strcmp(argv[k], "-v")==0 && k+1!=argc) {
			params.Verbosity = atoi(argv[++k]);
		} else if(::strcmp(argv[k], "-l")==0 && k+1!=argc) {
			leaveFile = argv[++k];
		}
	}
	return 0;
}

int OutputSettings(struct DH_TYPE::SPParameters& params, std::string& GTPath, int& task, char*& weightFile, char*& leaveFile) {
	std::cout << "Settings:" << std::endl;
   std::cout << "  (-t) Task:                 " << task << std::endl;
   std::cout << "  (-d) GTPath:               " << GTPath << std::endl;
   std::cout << "  (-w) WeightFile:           " << ((weightFile==NULL)?"":weightFile) << std::endl;
   std::cout << "  (-l) LeaveFile:            " << ((leaveFile==NULL)?"":leaveFile) << std::endl;
	std::cout << "  (-a) ArmijoIterations:     " << params.ArmijoIterations << std::endl;
	std::cout << "  (-b) BetheCountingNumbers: " << params.BetheCountingNumbers << std::endl;
	std::cout << "  (-c) C:                    " << params.C << std::endl;
	std::cout << "  (-i) CCCPIterations:       " << params.CCCPIterations << std::endl;
	std::cout << "  (-f) CRFEraseMessages:     " << params.CRFEraseMessages << std::endl;
	std::cout << "  (-j) CRFIterations:        " << params.CRFIterations << std::endl;
	std::cout << "  (-k) CRFMPIterations:      " << params.CRFMPIterations << std::endl;
	std::cout << "  (-x) CRFOuterExchange:     " << params.CRFOuterExchange << std::endl;
	std::cout << "  (-o) CRFOuterIterations:   " << params.CRFOuterIterations << std::endl;
	std::cout << "  (-e) epsilon:              " << params.epsilon << std::endl;
	std::cout << "  (-m) MPIterations:         " << params.MPIterations << std::endl;
	std::cout << "  (-p) p:                    " << params.p << std::endl;
	std::cout << "  (-y) ReadBinary:           " << params.ReadBinary << std::endl;
	std::cout << "  (-r) ReuseMessagesForF:    " << params.ReuseMessagesForF << std::endl;
	std::cout << "  (-v) Verbosity:            " << params.Verbosity << std::endl;
	return 0;
}

int main(int argc, char** argv) {
#ifdef USE_ON_WINDOWS
	std::string GTPath("..\\Data01");
#else
	std::string GTPath("Data01");
#endif
	std::string pattern("*.feature");

	struct DH_TYPE::SPParameters SPSettings;
	SPSettings.CRFIterations = 20;
	SPSettings.ArmijoIterations = 50;
	SPSettings.CRFEraseMessages = 1;
	SPSettings.ReuseMessagesForF = 0;
	SPSettings.argc = argc;
	SPSettings.argv = argv;
	SPSettings.epsilon = 1.0;
	SPSettings.CRFOuterIterations = 200;
	SPSettings.CRFMPIterations = 50;
	SPSettings.C = 1.0;
	SPSettings.p = 2.0;
	SPSettings.MPIterations = 50;
	SPSettings.CRFOuterExchange = size_t(-1);
	SPSettings.BetheCountingNumbers = 0;
	SPSettings.CCCPIterations = 1;
	SPSettings.ReadBinary = 0;

	int task = 0;
	char* weightFile = NULL;
	char* leaveFile = NULL;
	ParseInput(argc, argv, SPSettings, GTPath, task, weightFile, leaveFile);
	TraverseDirectory(GTPath, pattern, false, SPSettings.FeatureFiles);
   std::sort(SPSettings.FeatureFiles.begin(), SPSettings.FeatureFiles.end());
	
	try {
		StructuredPrediction<NUM_TYPE,SP_TYPE,DH_TYPE> SP(&SPSettings);

		if(SP.GetRank()==0) {
			OutputSettings( SPSettings, GTPath, task, weightFile, leaveFile );
		}

		std::vector<NUM_TYPE> w(SP.GetNumFeatures(), NUM_TYPE(1.0));
		if(weightFile!=NULL) {
			std::ifstream ifs(weightFile, std::ios_base::binary);
			int numFeatures, idTask;
			ifs.read((char*)&idTask, sizeof(int));
         if((idTask==0) || (idTask==1)) {
            ifs.read((char*)&numFeatures, sizeof(int));
         } else {
            numFeatures = idTask;
         }
			if(size_t(numFeatures)==w.size()) {
				std::vector<double> wf(numFeatures, 0.0);
				ifs.read((char*)&wf[0], numFeatures*sizeof(double));
				for(int k=0;k<numFeatures;++k) {
					w[k] = NUM_TYPE(wf[k]);
				}
			} else {
				std::cout << "Weight file dimension does not match." << std::endl;
			}
			ifs.close();
		}

      std::ofstream ofs;
      if(leaveFile!=NULL && SP.GetRank()==0) {
         ofs.open(leaveFile, std::ios_base::out | std::ios_base::binary);
         ofs.write((char*)&task, sizeof(int));
      }

		if(task<2) {
			SP.Iterate(&w);

			if(SP.GetRank()==0) {
				if(leaveFile!=NULL) {
					int numFeat = int(w.size());
					ofs.write((char*)&numFeat, sizeof(int));
				}
				std::cout << "w = [";
				for(size_t k=0;k<w.size();++k) {
					std::cout << w[k] << " ";
					if(leaveFile!=NULL) {
						double tmp = w[k];
						ofs.write((char*)&tmp, sizeof(double));
					}
				}
				std::cout << "]" << std::endl;
			}
		}

		if(task>0) {
			std::vector<std::map<size_t, std::vector<NUM_TYPE> > > Beliefs;
			SP.Predict(&w, Beliefs);

			if(SP.GetRank()==0) {
				if(leaveFile!=NULL) {
					int numSamples = int(Beliefs.size());
					ofs.write((char*)&numSamples, sizeof(int));
				}
				for(std::vector<std::map<size_t, std::vector<NUM_TYPE> > >::iterator x=Beliefs.begin(),x_e=Beliefs.end();x!=x_e;++x) {
					//std::cout << "Sample " << x-Beliefs.begin() << ":" << std::endl;
					if(leaveFile!=NULL) {
						int numRegions = int(x->size());
						ofs.write((char*)&numRegions, sizeof(int));
					}
					for(std::map<size_t, std::vector<NUM_TYPE> >::iterator r=x->begin(),r_e=x->end();r!=r_e;++r) {
						//std::cout << "  " << r->first << ":";
						if(leaveFile!=NULL) {
							int belSize = int(r->second.size());
							ofs.write((char*)&belSize, sizeof(int));
						}
						for(std::vector<NUM_TYPE>::iterator rx=r->second.begin(),rx_e=r->second.end();rx!=rx_e;++rx) {
							//std::cout << " " << *rx;
							if(leaveFile!=NULL) {
								double bel = double(*rx);
								ofs.write((char*)&bel, sizeof(double));
							}
						}
						//std::cout << std::endl;
					}
				}
			}
		}

		if(leaveFile!=NULL) {
			ofs.close();
		}

		SP.Clear();
	} catch(int e) {
		std::cout << "Error: " << e << std::endl;
	}

	return 0;
}
