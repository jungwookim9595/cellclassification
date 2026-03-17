/** 
    This code is a C wrapper, through which fortran code can invoke 
    ANN c++ code to get the wall distance.
*/
/**
   Author: Xudong Xiao, Ph.D. of Aerospace Engineering, NCSU
   Date : Nov-22-2003
*/


#include<fstream>
#include<iostream>
#include<ANN/ANN.h>
#include<stdlib.h>
#include<cmath>
static ANNpointArray  surf_pts;      // surface points

static ANNkd_tree   *the_tree;      // search structure
using namespace std;
#include<string.h>
extern "C" {

  /// initialize all local static variables of this module
  void wldst_init_(char *filename,int length){

    ifstream ifs;
    char *str=new char[length+1];
    strncpy(str,filename,length);
    str[length]=0;
    cout<<"File name="<<str<<"  "<<length<<endl<<flush;
    ifs.open(str);

    int num; // number of surface points
    ifs>>num;
    if(ifs.bad()){
      cerr<<"File "<<str<<" doesn't exist."<<endl<<flush;
      exit(0);
    }
   
    surf_pts = annAllocPts(num,3);    // 
    

    /// read surface points
    for(int i=0;i<num;i++){
      ifs>>surf_pts[i][0]
         >>surf_pts[i][1]
         >>surf_pts[i][2];
      if(!ifs.good()){
        cerr<<"Reading surface pts file failed."<<endl<<flush;
        exit(-1);
      }

    }

    ifs.close();

    // initialize the kd-tree with surface points
    the_tree = new ANNkd_tree(        // build search structure
                              surf_pts, // the data points
                              num,      // number of points
                              3);       // dimension of space
    delete str;
  }
 
  void wldst_create_tree_(double *x, double *y, double *z, const int *num)
  {
    surf_pts = annAllocPts(*num,3);    // 
    for(int i=0;i<*num;i++){
      surf_pts[i][0] = x[i];
      surf_pts[i][1] = y[i];
      surf_pts[i][2] = z[i];
    }
    the_tree = new ANNkd_tree(        // build search structure
                              surf_pts, // the data points
                              *num,      // number of points
                              3);       // dimension of space
    return;
  }
  
  void wldst_delete_tree_(){
    delete the_tree;
    the_tree = 0;
    annDeallocPts(surf_pts);
    surf_pts = 0;
  }
 
  void wldst_getdist_(ANNpoint query_pt, double *dist){

    static const int knear=5;      // the number of nearest points
    static const double eps = 0.;  // tolerance
    static ANNidx  nn_idx[knear];      /// near neighbor indices
    static ANNdist dist_array[knear];      // near neighbor distances

    the_tree->annkSearch(          // search
                         query_pt,     // query point
                         knear,        // number of near neighbors
                         nn_idx,       // nearest neighbors (returned)
                         dist_array,   // distance (returned)
                         eps);         // error bound
    *dist = sqrt(dist_array[0]);

    for(int i=0;i<knear;i++){
      cout<<"dist["<<i<<"]="<<sqrt(dist_array[i])<<" indx="<<nn_idx[i]<<endl;
    }
    return;
  }

  void wldst_getdist_xyz_(double *x,
                          double *y,
                          double *z,
                          int *indx, 
                          double *dist)
  {

    static const int knear=5;      // the number of nearest points
    static const double eps = 0.;  // tolerance
    static ANNidx  nn_idx[knear];      /// near neighbor indices
    static ANNdist dist_array[knear];      // near neighbor distances
    static double pt[3];
    pt[0] = *x;
    pt[1] = *y;
    pt[2] = *z;
    the_tree->annkSearch(          // search
                         pt,     // query point
                         knear,        // number of near neighbors
                         nn_idx,       // nearest neighbors (returned)
                         dist_array,   // distance (returned)
                         eps);         // error bound
    *dist = sqrt(dist_array[0]);
    *indx = nn_idx[0];
   
    return;
  }

  void wldst_getdist_xyzs_(double *x,
                           double *y,
			   double *z, double *normal,int *nsize,
                           int *indx, 
                           double *dist,double *Sign)
  {

    static const int knear=41;      // the number of nearest points
    static const double eps = 0.;  // tolerance
    static ANNidx  nn_idx[knear];      /// near neighbor indices
    static ANNdist dist_array[knear];      // near neighbor distances
    static double pt[3];
    int nn = *nsize;
    pt[0] = *x;
    pt[1] = *y;
    pt[2] = *z;
    the_tree->annkSearch(          // search
                         pt,     // query point
                         knear,        // number of near neighbors
                         nn_idx,       // nearest neighbors (returned)
                         dist_array,   // distance (returned)
                         eps);         // error bound
    *dist = sqrt(dist_array[0]);
    *indx = nn_idx[0];
    int ipos,ineg;  // two counters for interior, exterior,
    ipos = ineg=0;
    for(int i=0;i<knear;i++){
      int k = nn_idx[i];
      int offsetx = (*nsize)*3+k;
      int offsety = (*nsize)*4+k;
      int offsetz = (*nsize)*5+k;
      double dr[3];
      dr[0] = *x - surf_pts[k][0];
      dr[1] = *y - surf_pts[k][1];
      dr[2] = *z - surf_pts[k][2];

      double dotp = dr[0]*normal[offsetx] + dr[1]*normal[offsety]+dr[2]*normal[offsetz];
      if (dotp>0) ipos ++;
      else ineg++;
    }
    if (ipos >ineg) *Sign = 1.0;
    else *Sign = -1.0;
   
    return;
  }

  void wldst_getdists_(double *x, double *y, double *z, double *d,


                       int *ii, int *jj, int *kk){
    static const int knear=5;      // the number of nearest points
    static const double eps = 0.;  // tolerance
    static ANNidx  nn_idx[knear];      /// near neighbor indices
    static ANNdist dist_array[knear];      // near neighbor distances
    int i0 = *ii;
    int j0 = *jj;
    int k0 = *kk;

    // translate column-major 3d fortran array to
    // 1d C array
    for(int k=1;k<=k0;k++){
      for(int j=1;j<=j0;j++){
	for(int i=1;i<=i0;i++){
	  double point[3];
	  int offxyz = i+1+ (j+1)*(i0+4)+(k+1)*(i0+4)*(j0+4);
	  int offd   = i-1+ (j-1)*(i0)+(k-1)*(i0)*(j0);
	  point[0] = x[offxyz];
	  point[1] = y[offxyz];
	  point[2] = z[offxyz];
	  the_tree->annkSearch(// search
			       point,// query point
			       knear,     // number of near neighbors
			       nn_idx,// nearest neighbors (returned)
			       dist_array,// distance (returned)
			       eps);// error bound
	  d[offd] = sqrt(dist_array[0]);
	}
      }
    }
  }


  void wldst_getdists_xyz_(double *x,
                           double *y,
                           double *z,
                           int *indx,
                           double *dist,int *knear)
  {
    // comment this line out -xxiao
    //static const int knear=100;      // the number of nearest points
    static const double eps = 0.;  // tolerance

    // define ANNidx and ANNdist as pointers - xxiao
    ANNidx*  nn_idx=new ANNidx[*knear];      /// near neighbor indices
    ANNdist* dist_array=new ANNdist[*knear];      // near neighbor distances
    static double pt[3];

    pt[0] = *x;
    pt[1] = *y;
    pt[2] = *z;
    the_tree->annkSearch(          // search
                         pt,     // query point
                         *knear,        // number of near neighbors
                         nn_idx,       // nearest neighbors (returned)
                         dist_array,   // distance (returned)
                         eps);         // error bound

    for(int inear=0;inear<*knear;inear++){
    *(indx + inear) = nn_idx[inear];
    *(dist + inear) = sqrt(dist_array[inear]); 
    } 
    // newly added two lines by xxiao
    delete []nn_idx; delete [] dist_array;
    return;
  }
}
