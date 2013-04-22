//~ Copyright 2013 Luciano Spinello
//~ This file is part of LIBGeFP.
//~ 
//~ LIBGeFP is free software: you can redistribute it and/or modify
//~ it under the terms of the GNU General Public License as published by
//~ the Free Software Foundation, either version 3 of the License, or
//~ (at your option) any later version. 
//~ 
//~ LIBGeFP is distributed in the hope that it will be useful,
//~ but WITHOUT ANY WARRANTY; without even the implied warranty of
//~ MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//~ GNU General Public License for more details.
 //~ 
//~ You should have received a copy of the GNU General Public License
//~ along with LIBGeFP.  If not, see <http://www.gnu.org/licenses/>.


/*! \mainpage 
 *
 * \section need Why do I need Geometrical FLIRT Phrases (GFP)?
 * Let's say you have a large 2D scan dataset and you need to localize, or to check for loop closing.
 * Usually, a brute forcing ICP scan matching one vs all solves the problem. 
 * In case you need to do this process online, or repetitively, this may be computationally expensive. 
 * 
 * GFP comes to help in this situations. GFP makes use of a FLIRT bag-of-word representation, to create "phrases": combinations of words that 
 * are robust to noise and are able to encode complicate sequential patterns. In contrast to other approaches, GFP exploits the fact that a scan
 *  is a 1D manifold to create a rotation invariant scan representation.
 * This allows for reliable and fast retreival of 2D laser scans. Way more reliable 
 * than standard bag-of-words approaches.
 * 
 * In practice, given a scan represented by FLIRT words eg:
 * 
 * \verbatim scan example: 50 12 56 87 94 52 68 \endverbatim
 * 
 * LIBGeFP matches it with a dataset containing 
 * \verbatim scan1: 13 12 56 8 94 52 687
scan2: 4 58 66 45 33 6 1
scan3: 1 223 3
scan4: 45 33 6 1 12 56
...  \endverbatim
 * \section flirtconnection How can I use GFPs?
 * <b>You can use it as a library for standalone operations (see executables on Example section) or you can use it as a plugin of FLIRTlib for global one-shot localization.</b>
 *
 * The latest version of <a href="http://openslam.org/flirtlib.html">FLIRTLib</a> fully supports LIBGeFP, and it will download LIBGeFP automatically at compile time.
 * 
 * 
 * \section intro Why are GFPs useful?
 *
 * Place recognition, i.e., the problem of recognizing
 * if the robot is navigating in an already visited place, is a
 * fundamental problem in mobile robot navigation. Efﬁcient
 * solutions to this problem are relevant for effectively localizing
 * robots and for creating maps in real time. Relatively few methods
 * have been proposed to efﬁciently solve this problem in very large
 * environments using 2D range data. <b>In this paper, we introduce
 * geometrical FLIRT phrases (GFPs) as a novel retrieval method
 * for very efﬁcient and precise place recognition. GFPs perform
 * approximate 2D range data matching, have low computational
 * cost, can handle complicated partial matching patterns and are
 * robust to noise</b>. Experiments carried out with publicly available
 * datasets demonstrate that GFPs largely outperform state-of-theart approaches in 2D range-based place recognition in terms of
 * efﬁciency and recall. We obtain retrieval performances with more
 * than 85% recall at 99% precision in less than a second, even on
 * data sets obtained from several kilometer long runs.
 * 
 * see the <a href="http://www.informatik.uni-freiburg.de/~spinello/tipaldiICRA13.pdf">paper</a> 
 * \section install_sec Installation
 * Download the library. The library relies on cmake to generate the Makefiles.
 * 
 * Go to the \c root directory of your project and run
 * \verbatim $ mkdir build
$ cd build  
$ cmake .. 
$ make \endverbatim
 Binaries are generated in \c build/bin and \c build/lib
 * 
 * Library can be installed in your system by running
 * \verbatim $ make install \endverbatim
 * 
 * 
 * 
 * The software depends on the following external libraries
 * \li <em> Boost >= 1.4 (special_functions/binomial) </em>
 * \section ex_sec Examples
 * Example files can be found in \c bin/. For testing a FLIRT words dataset has been included in \c data_example/
 * 
 * 
 * 
 * \verbatim bin/gefp_cl \endverbatim reads a dataset where each scan is composed of FLIRT words and retreives best matches in the dataset for benchmarking. 
 * Many parameters can be selected from command line (kernel size, bag-of-words/distances, etc)
 * 
 * \verbatim bin/gefp_cl_onequery \endverbatim same as above but with one scan. This is an example for understanding how to use it.
 * 
 * \section ref How to cite LIBGeFP
 *  "Geometrical FLIRT Phrases for Large Scale Place Recognition in 2D Range Data", G. D. Tipaldi, L. Spinello, W. Burgard -- Int. Conf. Robotics and Automation (ICRA) 2013
 *  
 * \section license License
 * GFP LIBGeFP Copyright (c) 2013 Luciano Spinello, licensed under GPL ver. 2.0
 * 
 * GFP is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version. 
 * GFP is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License
 * along with GFP.  If not, see <http://www.gnu.org/licenses/>.
 */
 
#ifndef GFP_NGN_H
#define GFP_NGN_H

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <set>
#include <float.h>
#include <limits.h>
#include <stdlib.h>
#include <stdio.h>
#include <boost/math/special_functions/binomial.hpp>
#include <sys/time.h>  
#include <algorithm>


//~ Basic and bag of distances  defaults   
#define DEFAULT_BOWDST_START 0
#define	DEFAULT_BOWDST_INTERVAL 0.2
#define DEFAULT_BOWDST_END 15
#define DEFAULT_BOWSUBTYPE 0
#define DEFAULT_ALPHASMOOTH 0.4
#define DEFAULT_BAGDISTANCE 0
#define DEFAULT_CACHEBINOMIAL 10000

/**
 * Contains a 2D scan represented by FLIRT words identified by their index, their (TF-IDF) weights, their norm for GFP
 * 
 * @author Luciano Spinello
 */					 
 
class scan_bow
{
	public:
		std::vector <int> w, word_weight_unnormalized;
		std::vector <double> w_x, w_y, word_weight, tfidf_w;
		double sum_weight,  norm_wgv;
		scan_bow(int no)
		{
			w= std::vector <int> (no);
			w_x = std::vector <double> (no);
			w_y = std::vector <double> (no);
		}
};
 
 

/**
 * Caches the word orders in a scan for GFP indexing
 * 
 * @author Luciano Spinello
 */	
class tf_idf_db_ordercache
{
	public:
		std::vector<int> pos;
};


/**
 * Contains TF-IDF weight for a  FLIRT word
 * 
 * @author Luciano Spinello
 */	
class tf_idf_db
{
	public:
		//~ per doc
		std::vector <tf_idf_db_ordercache> word_order;
		std::vector <int> doc_id;
		std::vector <int> term_count_unnormalized;
		std::vector <double> tf_idf_doc_normed;
		std::vector <double> ntf_idf_doc_normed;
		std::vector <double> wf_idf_doc_normed;
		std::vector <int> num_words;
		std::vector <double> term_count;
		
		//~ per term
		int num_doc_containing_the_word, corpus_size;
		double idf;
		
		tf_idf_db()
		{
			num_doc_containing_the_word = 0;
			corpus_size = 0;
		}
};

/**
 * Geometrical FLIRT Phrases (GFP) for matching 2D laser scans represented FLIRT words
 * 
 * It includes methods for building the search index and matching
 * 
 * <a href="http://www.informatik.uni-freiburg.de/~spinello/tipaldiICRA13.pdf">"Geometrical FLIRT Phrases for Large Scale Place Recognition in 2D Range Data", G. D. Tipaldi, L. Spinello, W. Burgard -- Int. Conf. Robotics and Automation (ICRA) 2013</a> 
 * @author Luciano Spinello
 */		
class gefp_engine
{
	private:
		//~ vars
		std::vector <scan_bow> laserscan_bow;
		std::vector < std::pair <double, int> > scoreset;
		std::vector <tf_idf_db> tf_idf;
		std::string fileoutput_rootname;
		int dictionary_dimensions, start_l, stop_l, max_bow_len, wgv_kernel_size, bow_type, bow_subtype;
		double anglethres, bow_dst_start, bow_dst_interval, bow_dst_end, alpha_vss;
		uint number_of_scans, kbest;
		std::vector<double> cached_binomial_coeff, mtchgfp_rc_idf_sum, normgfp_rc_idf_sum;
		std::vector <int> mtchgfp_min_det_idx, mtchgfp_max_det_idx, mtchgfp_rc_weak_match, normgfp_rc_weak_match;
		std::vector<char> mtchgfp_used_doc_idx;

		//~ functions
 		double norm_gfp(std::vector <int> & query_v);
 		void matching_bow(std::vector <int> &query_v );
		void matching_gfp(std::vector <int> &query_v );
		void voting_tfidf_weak_verificationOLD(std::vector <int> &query_v );		
		void reformulate_to_bagofdistances(void);
		void cache_binomial_coeff(void);
	
	public:

		//~ functions		
		/**
		 * Reads file generated by FLIRTLIB in which each scan is described as a sequence of FLIRT words, represented each by a number
		 * 
		 * @author Luciano Spinello
		 */		 
 		int read_wordscan_file(std::string filename);

		/**
		 * Inserts a scan described as a sequence of FLIRT words, represented each by a number 
		 * @param wordscan scan identified as a sequence of ids
		 * @param xpos,ypos metric position of each word in \c wordscan
		 * @author Luciano Spinello
		 */		 
 		void insert_wordscan(std::vector <int> wordscan, std::vector <double> xpos, std::vector <double> ypos);
 		

		/**
		 * Builds TF-IDF index for standard and weak verification matching methods
		 * 
		 * It implements IDF, TF, TF-IDF, wordcount and improved TF-IDF models proposed in:
		 * <a href="http://comminfo.rutgers.edu/~muresan/IR/Docs/Articles/ipmSalton1988.pdf"> Gerard Salton and Christopher Buckley: "Term-weighting approaches in automatic text retrieval", Information processing & management, vol. 24, no. 5, 1988, Elsevier </a>
		 * @author Luciano Spinello
		 */		
 		void build_tfidf(void);

		/**
		 * Matches all the scans in the dataset vs all the scans in the dataset
		 * 
		 * Saves the k-best results on disk for each query along with computational time 
		 * 
		 * @param dtype  kind of matching method: 1 standard bag-of-words, 2 geometrical FLIRT phrases
		 * @author Luciano Spinello
		 */
		void run_evaluation(int dtype);
		

		/**
		 * Matches a query scan with the dataset
		 * 
		 * Computes the k-best results on disk
		 * 
		 * @param dtype  kind of matching method: 1 standard bag-of-words, 2 geometrical FLIRT phrases
		 * @param query_v a query scan, composed by a vector of numbers, each indicating a FLIRT word
		 * @param scoreoutput a pointer to a sorted vector of pairs containing <scorematch, index of the scan in the dataset>. 
		 * @author Luciano Spinello
		 */
		void query(int dtype, std::vector <int>   &query_v, std::vector < std::pair <double, int> > **scoreoutput);


		/**
		 * Prepares indeces and cache for matching. Executed once at the beginning.
		 * Builds TF-IDF index for the dataset and norms all vectors on the dataset. Allocates also needed memory.
		 * Optionally it generates bag-of-distances
		 * @author Luciano Spinello
		 */
		void prepare(void);

		/**
		 * Constructor
		 * 
		 * @param krnl Kernel size
		 * @param kbt # of k-best results
		 * @param bt 1 for bag-of-distances, 0 otherwise
		 * @param bstype flavor of TF-IDF in case of standard bag-of-words: 0 standard TFIDF, 1 sublinear TFIDF scaling, 2 lenght smoothing TFIDF, see \link gefp_engine::build_tfidf\endlink
		 * @param a_vss alpha_smoothing in case of standard bag-of-words with lenght smoothing TFIDF (0.4 default)
		 * @author Luciano Spinello
		 */  

		gefp_engine (int krnl, int kbt, int bt=DEFAULT_BAGDISTANCE, int bstype=DEFAULT_BOWSUBTYPE, double a_vss=DEFAULT_ALPHASMOOTH)
		{
			bow_type = bt;
			kbest = kbt;
			wgv_kernel_size = krnl;
 			bow_subtype= bstype;
			alpha_vss = a_vss;
			
			//~ basic defaults for bag of distances
			bow_dst_start= DEFAULT_BOWDST_START;
		    bow_dst_interval=DEFAULT_BOWDST_INTERVAL;
		    bow_dst_end=DEFAULT_BOWDST_END;

		}
		
};

/**
 * Simple string tokenizer (avoids additional dependencies)
 * 
 * @author Luciano Spinello
 */	
void LSL_stringtoken(const std::string& str, std::vector<std::string>& tokens, const std::string& delimiters);
bool isBettermatched(std::pair <double, int> x, std::pair <double, int> y); 
#endif
