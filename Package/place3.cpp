/*
 * Copyright 2012-2013 TopCoder, Inc.

 *
 * This code was developed under U.S. government contract NNH10CD71C. 

 *
 * Licensed under the Apache License, Version 2.0 (the "License");

 * You may not use this file except in compliance with the License.

 * You may obtain a copy of the License at:
 *     http://www.apache.org/licenses/LICENSE-2.0

 * Unless required by applicable law or agreed to in writing, software

 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.

 * See the License for the specific language governing permissions and

 * limitations under the License.
 */

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
#include <sstream>
#include <cmath>
#include <sys/time.h>
#include <utility>
#include <map>
#include <stack>
#include <cstring>
#include <cstdint>

using namespace std;
#define W 1236
#define H 630
#define R 40
#define R2 1600  
//~ #define VERBOSE
//~ #define PRINTRAM

#define MAX(x,y) ((x>y)?(x):(y))
#define MIN(x,y) ((x<y)?(x):(y))
typedef vector< vector<int> > vvi;
typedef vector<int> vi;
typedef vector< vector<double> > vvd;
typedef vector<double> vd;
typedef vector< vector<bool> > vvb;
typedef vector<bool> vb;
typedef vector< vector<string> > vvs;
typedef vector<string> vs;
typedef vector< vector<uint8_t> > vvui8;
typedef vector<uint8_t> vui8;
typedef map<string,int>  mapi;
typedef pair<int,int> pii;

static double t0;
static bool skip_training = false;
static bool modern;

class Grid{
  public: 
  int x0, y0, dx, dy, nx, ny;
  vector<int> x;
  vector<int> y; 
  int txpos(int xpos){ return MIN( nx-1, MAX(0, (xpos-x0)/dx ) ); }
  int typos(int ypos){ return MIN( ny-1, MAX(0, (ypos-y0)/dy ) ); }
  void set(int _x0, int _y0, int _dx, int _dy)
  {
    x0 = _x0;
    y0 = _y0;
    dx = _dx;
    dy = _dy;
    nx = (W - 2*x0)/dx +2;
    ny = (H - 2*y0)/dy +2;
    // center grid:
    x0 = (W - (nx-1)*dx)/2;
    y0 = (H - (ny-1)*dy)/2;
    x.resize(nx);
    for(int i=0; i<nx; ++i)
      x[i] = x0 + i*dx;
    y.resize(ny);
    for(int i=0; i<ny; ++i)
      y[i] = y0 + i*dy;       
  }
};

class Frame {
  public:
  int y0;
  int y1;
  int x0;
  int x1;
  Frame():y0(0),y1(0),x0(0),x1(0){};
  Frame(int _y0, int _y1, int _x0, int _x1):y0(_y0),y1(_y1),x0(_x0),x1(_x1){};
  void init(int _y0, int _y1, int _x0, int _x1){ 
    y0 = _y0; y1 = _y1; x0 = _x0; x1 = _x1;
  }
  void init(Frame& f) { 
    y0 = f.y0; y1 = f.y1; x0 = f.x0; x1 = f.x1;
  }
  void update(int y, int x) {
    y0 = MIN(y0, y);
    y1 = MAX(y1, y);
    x0 = MIN(x0, x);
    x1 = MAX(x1, x);
  }
  int get_height() { return (y1-y0); }
  int get_width() { return (x1-x0); }
};

///----------
/// UTILTY
///

double getTime()
{
    timeval tv;
    gettimeofday(&tv, 0);
    return tv.tv_sec + tv.tv_usec * 1e-6;
}
void  initTime(){ t0 = getTime(); }
double Time(){ return getTime()-t0; }

#ifdef PRINTRAM
inline void printRam()
{
 struct rusage ru;
 getrusage(RUSAGE_SELF, &ru);
 cerr << "ram (Mbytes): " << ru.ru_maxrss/1000 << endl;
}
#else
inline void printRam() {}
#endif

template<class T>
void clear_vv(vector< vector<T> >& t)
{
  if(!t.empty()) {
    for(int i=0, isize=t.size(); i<isize; ++i)
      vector<T>().swap( t[i] );
    vector< vector<T> >().swap( t );
  }
}

string int2string(int i) {
  string ret = "";
  if(i>=1000) {
    ret += ( (i/1000) + '0' );
    i %= 1000;
  }
  if(i>=100) {
    ret += ( (i/100) + '0' );
    i %= 100;
  }
  if(i>=10) {
    ret += ( (i/10) + '0' );
    i %= 10;
  }    
  ret += ( i + '0' );
  return ret;
}

int distance2(int qy, int qx, int py, int px) {
    return ( (qy-py)*(qy-py)+(qx-px)*(qx-px) );
  }

void browse(int i, int j, vvb& visited, vvb& condition, Frame& f, int& counter)
{
    stack<pii> s;
    s.push(pii(i,j));
    while(!s.empty()) {
      pii pixel = s.top();
      ++counter;
      int ii = pixel.first;
      int jj = pixel.second;
      visited[ii][jj] = true;
      f.update(ii, jj);
      s.pop();
      if( ii<H-1 && condition[ii+1][jj] && !visited[ii+1][jj] )
	s.push(pii(ii+1,jj));
      if( ii>0 && condition[ii-1][jj] && !visited[ii-1][jj] )
	s.push(pii(ii-1,jj));
      if( jj<W-1 && condition[ii][jj+1] && !visited[ii][jj+1] )
	s.push(pii(ii,jj+1));
      if( jj>0 && condition[ii][jj-1] && !visited[ii][jj-1] )
	s.push(pii(ii,jj-1));
    }
}

//----------------------------------------------------------------------
// GLOBAL FUNCTIONS

template<class T>
void extrapolate_border(vector< vector<T> >& mat)
{
  for(int i=1; i<H-1; ++i) {
    mat[i][0] = mat[i][1];
    mat[i][W-1] = mat[i][W-2];
  }
  for(int j=1; j<W-1; ++j) {
    mat[0][j] = mat[1][j];
    mat[H-1][j] = mat[H-2][j];
  }
  mat[0][0] = (mat[0][1] + mat[1][0])/2;	
  mat[H-1][0] = (mat[H-2][0] + mat[H-1][1])/2;	
  mat[0][W-1] = (mat[1][W-1] + mat[0][W-2])/2;	
  mat[H-1][W-1] = (mat[H-2][W-1] + mat[H-1][W-2])/2;
}

void cumulative(vvd& original, vvd& cumul) {
    cumul.assign(H, vd(W, 0.0) );
    for(int j=0; j<W; ++j)
      cumul[0][j] = original[0][j];
    for(int i=1; i<H; ++i)
      for(int j=0; j<W; ++j)
	cumul[i][j] = cumul[i-1][j] + original[i][j];
    for(int j=1; j<W; ++j)
      for(int i=0; i<H; ++i)
	cumul[i][j] += cumul[i][j-1];
}

void cumulative2(vvd& original, vvd& cumul) {
    cumul.assign(H, vd(W, 0.0) );
    for(int j=0; j<W; ++j)
      cumul[0][j] = original[0][j]*original[0][j];
    for(int i=1; i<H; ++i)
      for(int j=0; j<W; ++j)
	cumul[i][j] = cumul[i-1][j] + original[i][j]*original[i][j];
    for(int j=1; j<W; ++j)
      for(int i=0; i<H; ++i)
	cumul[i][j] += cumul[i][j-1];
}

double cellsum(vvd& cumul, Frame& f) {
  // y0/x0 index included , y1/x1 index excluded
  int y0=f.y0, y1=f.y1, x0=f.x0, x1=f.x1;
  if(x0==x1 || y0==y1)
    return 0.0;
  if(x0==0 || y0==0) {
    if(x0==0 && y0==0)   
      return cumul[y1][x1];
    if(x0==0) {
      --x1, --y0, --y1;
      return (cumul[y1][x1] -cumul[y0][x1]);
    }
    if(y0==0) {
      --x0, --x1, --y1;
      return (cumul[y1][x1] -cumul[y1][x0]);
    }
  }
  --x0, --x1, --y0, --y1;
  return (cumul[y0][x0] + cumul[y1][x1] - cumul[y1][x0] -cumul[y0][x1]);
}

// END OF GLOBAL FUNCTION
//----------------------------------------------------------------------

class Point
{
  public:
  int y;
  int x;
  
  Point():y(0),x(0){}
  Point(int _y, int _x):y(_y), x(_x){}
  int distance2(int py, int px) {
    return ( (y-py)*(y-py)+(x-px)*(x-px) );
  }
};


class Annotator
{
    public:
    int nbefore;
    int nlifetime;
    int time;
    vector<Point> target;
    vector<Point> other;
    vvb tstate;
    
    Annotator():nbefore(0), nlifetime(0), time(0){};
    Annotator(int nb, int nl, int ti):nbefore(nb), nlifetime(nl), time(ti){};
    void build_tstate(Grid& grid)
    {
      tstate.assign(grid.ny, vb(grid.nx, false) );
      for(auto s=target.begin(), send=target.end(); s!=send; ++s) {
	int i0 = grid.txpos(s->x-R), i1 = grid.txpos(s->x+R);
	int j0 = grid.typos(s->y-R), j1 = grid.typos(s->y+R);
	for(int j=j0; j<=j1; ++j) {
	  int y = grid.y0 + j*grid.dy;	
	  for(int i=i0; i<=i1; ++i) {
	    int x = grid.x0 + i*grid.dx;
	    if(!tstate[j][i] && s->distance2(y, x)<=R2) 
	      tstate[j][i] = true;
	  }
	}
      }
    }
    
    void contrib(Grid& grid, vvd& tconf)
    {
      for(int j=0; j<grid.ny; ++j)	
	for(int i=0; i<grid.nx; ++i)
	  if(tstate[j][i]) tconf[j][i] += 1.0;
    }    

  bool on_target(int y, int x) 
  {
    for(auto s=target.begin(), send=target.end(); s!=send; ++s)
      if( s->distance2(y,x)<=R2 ) return true;
    return false;
  }
  bool on_other(int y, int x) 
  {
    for(auto s=other.begin(), send=other.end(); s!=send; ++s)
      if( s->distance2(y,x)<=R2 ) return true;
    return false;
  }
    
};


class Special : public Point
{
  public:
  int type;
  int size;
  bool active;  
  double conf;
  double score;
  Frame f;
  Special():Point(){}
  Special(int _y, int _x, int _type, int _size):
    Point(_y, _x), type(_type), size(_size), active(true){}
};

class Label : public Point
{
  public:
  double proba;
  Label(int _y, int _x, double _proba): Point(_y, _x), proba(_proba){}  
};

bool high_score(const Special& a, const Special& b) {
  return (a.score>b.score);
}

bool high_proba(const Label& a, const Label& b) {
  return (a.proba>b.proba);
}

bool high_pii_second(const pii& a, const pii& b) {
  return (a.second>b.second);
}

class BadBlob 
{
  public:
  Frame f;
  BadBlob() {}
  BadBlob(Frame & _f) { f.init(_f); }
  void init(Frame &_f){ f.init(_f); }
  bool intersect( int y, int x, int shift) {
    if( y<f.y0-shift || y>f.y1+shift || x<f.x0-shift || x>f.x1+shift )
      return false;
    return true;
  }
};


class Image
{
  public:
    string id;
    int status;
    vector<Annotator> annotator;
    int color_min;
    int color_max;
    int vote_modern;
    int vote_other;
    bool call_for_annotations;
    int processed_annotator;
    
    // 1236*630
    vvui8 rdata;
    vvui8 gdata;
    vvui8 bdata;
    vvd grayi;
    vvd gray;
    vvd grad;
    vvd harris;
    vvd cumgrad;
    vvd cumgrad2;
    vvd hue;
    vvd value;
    vvd saturation;
    vector<Special> special;
    vector<Label> label;
    vector<BadBlob> badb;
    // tiled value
    //~ vvd tgradi;
    //~ vvd tgradv;
    
    Image():id(""){};
    
  void receive_data(vi& imgdata)
  {
    #ifdef VERBOSE
      cerr << "start receive_data" << endl;
    #endif
    // global init
    badb.clear();
    processed_annotator = 0;
    call_for_annotations = false;
    vote_modern = 0;
    vote_other = 0;
    label.clear();
    // rgb
    rdata.assign(H, vui8(W, 0) );  
    gdata.assign(H, vui8(W, 0) );  
    bdata.assign(H, vui8(W, 0) );      
    int iter = 2;
    color_min=255, color_max=0;
    for(int i=0; i<H; ++i)
      for(int j=0; j<W; ++j) {
	int d = imgdata[iter++];
	uint8_t c = (d & 0x00FF0000)>>16;
	rdata[i][j] = c;
	color_min = MIN(color_min, c), color_max=MAX(color_max, c);
	c = (d & 0x0000FF00)>>8;
	gdata[i][j] = c;
	color_min = MIN(color_min, c), color_max=MAX(color_max, c);
	c = (d & 0x000000FF);     
	bdata[i][j] = c;
	color_min = MIN(color_min, c), color_max=MAX(color_max, c);
      }
  }
  
  void normalize_color()
  {
    #ifdef VERBOSE
      cerr << "start normalize_color" << endl;
    #endif
    if(color_min!=color_max){
	uint8_t delta = color_max - color_min;
	for(int i=0; i<H; ++i)
	  for(int j=0; j<W; ++j) {
	    rdata[i][j] = (255*(rdata[i][j]-color_min))/delta;;
	    gdata[i][j] = (255*(gdata[i][j]-color_min))/delta;;
	    bdata[i][j] = (255*(bdata[i][j]-color_min))/delta;;
	}
    }
  }

    void build_cell(Grid& grid, int c, vvd& _gradi, vvd& _gradv) {
      #ifdef VERBOSE
	cerr << "start build_cell" << endl;
      #endif
      // Size tiles (crop to image's limit)
      int cell_hdx = c*grid.dx/2;
      int cell_hdy = c*grid.dy/2;
      vector<int> y0(grid.ny);
      vector<int> y1(grid.ny);
      vector<int> x0(grid.nx);
      vector<int> x1(grid.nx);
      for(int i=0; i<grid.ny; ++i) {
	y0[i] = MIN(H, MAX(0, grid.y[i] -cell_hdy) );
	y1[i] = MIN(H, MAX(0, grid.y[i] +cell_hdy ) );	
      }
      for(int j=0; j<grid.nx; ++j) {
	x0[j] = MIN(W, MAX(0, grid.x[j] -cell_hdx ) );
	x1[j] = MIN(W, MAX(0, grid.x[j] +cell_hdx ) );      
      }
      // assign grid point with avg over tile
      _gradi.assign(grid.ny, vd(grid.nx, 0.0));
      _gradv.assign(grid.ny, vd(grid.nx, 0.0));
      for(int i=0; i<grid.ny; ++i)
	for(int j=0; j<grid.nx; ++j) {
	  double factor = (double)((x1[j]-x0[j])*(y1[i]-y0[i]));
	  Frame f(y0[i], y1[i], x0[j], x1[j]);
	  if(factor!=0.0) {
	    factor = 1.0/factor;
	    _gradi[i][j] = factor*cellsum(cumgrad, f);
	    _gradv[i][j] = factor*cellsum(cumgrad2, f);
	    _gradv[i][j] -= _gradi[i][j]*_gradi[i][j];
	  }
	}
    }
    
    void build_harris_max(Grid& grid, int c, vvd& _hmax) {
      #ifdef VERBOSE
	cerr << "start build_harris_max" << endl;
      #endif
      // Size tiles (crop to image's limit)
      int cell_hdx = c*grid.dx/2;
      int cell_hdy = c*grid.dy/2;
      vector<int> y0(grid.ny);
      vector<int> y1(grid.ny);
      vector<int> x0(grid.nx);
      vector<int> x1(grid.nx);
      for(int i=0; i<grid.ny; ++i) {
	y0[i] = MIN(H, MAX(0, grid.y[i] -cell_hdy) );
	y1[i] = MIN(H, MAX(0, grid.y[i] +cell_hdy ) );	
      }
      for(int j=0; j<grid.nx; ++j) {
	x0[j] = MIN(W, MAX(0, grid.x[j] -cell_hdx ) );
	x1[j] = MIN(W, MAX(0, grid.x[j] +cell_hdx ) );      
      }
      // assign grid point with avg over tile
      _hmax.assign(grid.ny, vd(grid.nx, 0.0));
      for(int i=0; i<grid.ny; ++i)
	for(int j=0; j<grid.nx; ++j) {
	  double& ref = _hmax[i][j];
	  ref = 0.0;
	  for(int ii=y0[i]; ii<y1[i]; ++ii)
	    for(int jj=x0[j]; jj<x1[j]; ++jj) {
	      ref = MAX(ref, harris[ii][jj]);
	}
    }
  }

    void build_harris_avg(Grid& grid, int c, vvd& _havg) {
      #ifdef VERBOSE
	cerr << "start build_harris_avg" << endl;
      #endif
      // Size tiles (crop to image's limit)
      int cell_hdx = c*grid.dx/2;
      int cell_hdy = c*grid.dy/2;
      vector<int> y0(grid.ny);
      vector<int> y1(grid.ny);
      vector<int> x0(grid.nx);
      vector<int> x1(grid.nx);
      for(int i=0; i<grid.ny; ++i) {
	y0[i] = MIN(H, MAX(0, grid.y[i] -cell_hdy) );
	y1[i] = MIN(H, MAX(0, grid.y[i] +cell_hdy ) );	
      }
      for(int j=0; j<grid.nx; ++j) {
	x0[j] = MIN(W, MAX(0, grid.x[j] -cell_hdx ) );
	x1[j] = MIN(W, MAX(0, grid.x[j] +cell_hdx ) );      
      }
      // assign grid point with avg over tile
      _havg.assign(grid.ny, vd(grid.nx, 0.0));
      for(int i=0; i<grid.ny; ++i)
	for(int j=0; j<grid.nx; ++j) {
	  double factor = (double)((x1[j]-x0[j])*(y1[i]-y0[i]));
	  Frame f(y0[i], y1[i], x0[j], x1[j]);
	  if(factor!=0.0) {
	    factor = 1.0/factor;
	    _havg[i][j] = factor*cellsum(harris, f);
	  }
	}
  }


    void score_special(Grid& g, vvb& tactive) {
      int counter = 0;
      // Restrict world of special points and "score"      
      for(auto s=special.begin(), send=special.end(); s<send; ++s) {
	s->score = s->size;
	for(int i=0, iend=badb.size(); i<iend; ++i)
	  if( s->active && badb[i].intersect(s->y, s->x, 20) )
	    s->active = false;
      }
      
      // Sort and place special points
      sort(special.begin(), special.end(), high_score);
      for(int i=0, isize=special.size(); i<isize && counter<3; ++i) {
	if(special[i].active) {
	  int sy = special[i].y;
	  int sx = special[i].x;
	  double proba = 0.25 + 0.25*(special[i].type)-0.1*counter;
	  if(modern)
	    label.push_back( Label(sy, sx, proba) );
	  //~ update_score_spe( label.back() );
	  ++counter;
	  special[i].active = false;
	  // check other specials: disable if intersect
	  for(int j=i+1, jsize=special.size(); j<jsize; ++j)
	    if( special[j].active && special[j].distance2(sy, sx)<=4*R2 )
	      special[j].active = false;
	  // check grid: disable if intersect
	  int j0 = g.typos(sy-2*R);
	  int j1 = MIN(g.ny, g.typos(sy+2*R) + 1);
	  int k0 = g.txpos(sx-2*R);
	  int k1 = MIN(g.nx, g.txpos(sx+2*R) + 1);
	  for(int j=j0; j<j1; ++j) 
	    for(int k=k0; k<k1; ++k)
	      if( tactive[j][k] && special[i].distance2(g.y[j], g.x[k])<=4*R2 )
		tactive[j][k] = false;
	}
      }
    }
    
    void score_tiles(Grid& g, vvb& tactive) {
      
      int scounter = label.size();
      // Restrict world of grid points and "score"
      
      vvd hmax;
      build_harris_max(g, 1, hmax);
      
      vector<Label> cc;
      for(int j=0; j<g.ny; ++j) 
	for(int k=0; k<g.nx; ++k) {
	  for(int i=0, iend=badb.size(); i<iend; ++i) {
	    if(tactive[j][k] && badb[i].intersect(g.y[j], g.x[k], 40) )
	      tactive[j][k] = false;
	  }
	  // build a pool of candidate	  
	  if(tactive[j][k]) {
	    if(modern)
	      cc.push_back(Label(j, k, hmax[j][k] ) );
	    else if( hmax[j][k]<1.5 && hmax[j][k]>0.0 )
	      cc.push_back(Label(j, k, hmax[j][k] ) );
	  }
	}

      // sort and place grid points
      sort(cc.begin(), cc.end(), high_proba);
      // reset since bad point already disabled
      tactive.assign(g.ny, vb(g.nx, true) ) ;
      for(int tcounter=0, i=0; tcounter<12-scounter && i<cc.size(); ++i) {
	if(tactive[cc[i].y][cc[i].x]) {
	  int ty = g.y[cc[i].y];
	  int tx = g.x[cc[i].x];
	  if(modern)
	    label.push_back( Label(ty , tx, MAX(0.0, 0.035-0.003*tcounter) ) );
	  else
	    //~ label.push_back( Label(ty , tx, MAX(0.0, 0.025-0.0015*tcounter) ) );
	    label.push_back( Label(ty , tx, MAX(0.0, 0.025-0.0015*tcounter) ) );
	  ++tcounter;
	  // update grid points
	  for(int j=0; j<g.ny; ++j) 
	    for(int k=0; k<g.nx; ++k)
	      if( tactive[j][k] && 
		  distance2( ty, tx, g.y[j], g.x[k])<=4*R2 )
		    tactive[j][k] = false;
	}
      }
    }
    
    void score_training(Grid& g) {
      #ifdef VERBOSE
	cerr << "start score_training" << endl;
      #endif
      vvb tactive;
      tactive.assign(g.ny, vb(g.nx, true) ) ;
      score_special(g, tactive);
      // received annotions to check our findings
      if( label.size() )
	call_for_annotations = true;
    }
    
    void score_live(Grid& g) {
      #ifdef VERBOSE
	cerr << "start score_live" << endl;
      #endif
      vvb tactive;
      tactive.assign(g.ny, vb(g.nx, true) ) ;
      score_special(g, tactive);
      score_tiles(g, tactive);
    }

    void update_score_spe() {
      int asize = annotator.size();
      if(asize!=0) {
	for(auto s=label.begin(), send=label.end(); s!=send; ++s) {
	  int y = s->y;
	  int x = s->x;	
	  int scount = 0, ocount = 0;
	  for(auto s=annotator.begin(), send=annotator.end(); s!=send; ++s) {
	    scount += s->on_target(y, x);
	    ocount += s->on_other(y, x);
	  }
	  vote_modern += scount;
	  vote_other += ocount;
	}
      }
    }

    void sinit() {
      // feature
      special.clear();
      find_whitespot();
      find_colorblob();
      find_darkblob();
      find_whiteblob();
      find_water();
    }
    
    void build_gray()
    {
      #ifdef VERBOSE
	cerr << "start build_gray" << endl;
      #endif
      gray.assign(H, vd(W, 0.0));
      grayi.assign(H, vd(W, 0.0));
      double lmin = 255.0;
      double lmax = 0.0;      
      for(int i=0; i<H; ++i)
	for(int j=0; j<W; ++j) {
	  int rr = rdata[i][j];
	  int gg = gdata[i][j];
	  int bb = bdata[i][j];
	  grayi[i][j] = gray[i][j] = 0.30*(double)rr + 0.59*(double)gg + 0.11*(double)bb;
	  lmax = MAX(lmax, gray[i][j]);
	  lmin = MIN(lmin, gray[i][j]);	  
	}
      // Normalize
      double factor = lmax - lmin;
      if(factor>=1.0) {
	factor = 255.0 / factor;
	for(int i=0; i<H; ++i)
	  for(int j=0; j<W; ++j)
	    gray[i][j] = factor*(gray[i][j]-lmin);
       } else {
	return;
      }
    }
    
    void build_grad()
    {
      #ifdef VERBOSE
	cerr << "start build_grad" << endl;
      #endif
      grad.assign(H, vd(W, 0.0));

      // harris: inspired by opencv
      vvd dydy(H, vd(W, 0.0));
      vvd dydx(H, vd(W, 0.0));            
      vvd dxdx(H, vd(W, 0.0));
      
      {
	double factor = 1.0/255.0;
	double factor2 = factor*factor;
	double grad0, grad1;	
	for(int i=1; i<H-1; ++i)
	  for(int j=1; j<W-1; ++j) {
	    // sobel aperture 3 (scale 4)
	    // grad0 = y
	    grad0 = 2*gray[i-1][j] +  gray[i-1][j-1] + gray[i-1][j+1] -
	      (2*gray[i+1][j] +  gray[i+1][j-1] +gray[i+1][j+1]);
	    // grad1 = x  
	    grad1 = 2*gray[i][j-1] +  gray[i-1][j-1] +gray[i+1][j-1] -
	      (2*gray[i][j+1] +  gray[i-1][j+1] +gray[i+1][j+1]);
	    grad[i][j] = sqrt(grad0*grad0 + grad1*grad1)*factor;
	    dydy[i][j] = grad0*grad0*factor2;
	    dydx[i][j] = grad0*grad1*factor2;
	    dxdx[i][j] = grad1*grad1*factor2;
	  }
      }
      extrapolate_border(grad);
      	
      extrapolate_border(dydy);	
      extrapolate_border(dydx);	
      extrapolate_border(dxdx);	
      
      // Harris corner
      vvd cumuldydy;
      vvd cumuldydx;
      vvd cumuldxdx;
      cumulative(dydy, cumuldydy);
      cumulative(dydx, cumuldydx);
      cumulative(dxdx, cumuldxdx);
      harris.assign(H, vd(W, 0.0));
      int shift = 4;
      double kharris = 0.04;
      for(int i=0; i<H; ++i)
	for(int j=0; j<W; ++j) {
	    Frame f( MAX(0, i-shift), MIN(H, i+1+shift),
	      MAX(0, j-shift), MIN(W, j+1+shift) );
	    double factor = (f.get_height()-1)*(f.get_width()-1);
	    factor = 1/factor;
	    double a = cellsum(cumuldydy,f)*factor;
	    double b = cellsum(cumuldydx,f)*factor;
	    double c = cellsum(cumuldxdx,f)*factor;
	    harris[i][j] = a*c - b*b - kharris*(a + c)*(a + c);
	}
      
      cumulative(grad, cumgrad);
      cumulative2(grad, cumgrad2);
    }
    
    
    
    void find_colorblob()
    {

      vvb colored(H, vb(W, false));
      vvb visited(H, vb(W, false));
   
      for(int i=0; i<H; ++i)
	for(int j=0; j<W; ++j)
	  if(hue[i][j]>=100 && hue[i][j]<=300 && saturation[i][j]>0.5 && value[i][j]>0.7)
	    colored[i][j] = true;
      
      // fill flood to get candidate
      for(int i=0; i<H; ++i)
	for(int j=0; j<W; ++j)
	  if(colored[i][j] && !visited[i][j]) {
	    Frame f(i,i,j,j);
	    int counter = 0;
	    browse(i, j, visited, colored, f, counter);
	    int h = f.get_height();
	    int w = f.get_width();
	    //~ if(	counter >=8 && counter<=3000 ) {
	    if(	counter >=25 && counter<=3000 ) {
	      int cy = f.y0 + h/2;
	      int cx = f.x0 + w/2;
	      double gradi, gradv;
	      {
		int shift=6;
		Frame fs( MAX(0, f.y0-shift), 
			  MIN(H, f.y1+1+shift), 
			  MAX(0, f.x0-shift), 
			  MIN(W, f.x1+1+shift) );
		int hs = fs.get_height();
		int ws = fs.get_width();
		gradi = cellsum(cumgrad,fs)/(double)(ws*hs);
		gradv = cellsum(cumgrad2,fs)/(double)(ws*hs);
		gradv -= gradi*gradi;
	      }	      
	      cy = MAX(30, MIN(H-31, cy )); 
	      cx = MAX(30, MIN(W-31, cx )); 
	      Special spe(cy, cx, 1, h*w);
	      spe.f.init(f);
	      special.push_back(spe);
	  }
	}
      
      
    }
    


    void find_water()
    {
      #ifdef VERBOSE
	cerr << "start find_water" << endl;
      #endif
      vvb white(H, vb(W, false));
      vvb visited(H, vb(W, false));

      //~ thresh00 =  245;
      for(int i=0; i<H; ++i)
	for(int j=0; j<W; ++j) {
	  //~ if(grayi[i][j]>=thresh00 && saturation[i][j]<=0.05) {
	  //~ if(grayi[i][j]>=thresh00) {
	  if(abs(hue[i][j]-170)<=20 && 
	      saturation[i][j]>=0.1 && 
	      (saturation[i][j]+value[i][j])>=0.9 && gray[i][j]>=150.0) {
	    white[i][j] = true;
	  }
	}

      // fill flood to get candidate
      for(int i=0; i<H; ++i)
	for(int j=0; j<W; ++j)
	  if(white[i][j] && !visited[i][j]) {
	    Frame f(i,i,j,j);
	    int counter = 0;
	    browse(i, j, visited, white, f, counter);
	    //~ int h = f.get_height();
	    //~ int w = f.get_width();
	    if(	counter >= 100 ) {
	      //~ cerr << "WATER C : " << f.y0 + h/2 << " " << f.x0 + w/2 << " " << h << " " << w << endl;
	      badb.push_back(BadBlob(f));
	  }
	}
    }
    
    void find_whitespot()
    {
      #ifdef VERBOSE
	cerr << "start find_whitespot" << endl;
      #endif
      vvb white(H, vb(W, false));
      vvb visited(H, vb(W, false));
      
      //~ // Threshold heuristic
      double thresh00 = 0.0;     
      for(int i=0; i<H; ++i)
	for(int j=0; j<W; ++j)
	    thresh00 += gray[i][j];
      thresh00 /= W*H;
      thresh00 =  ( thresh00 + 0.5*(255.0-thresh00) );
      //~ thresh00 =  245;
      for(int i=0; i<H; ++i)
	for(int j=0; j<W; ++j) {
	  //~ if(grayi[i][j]>=thresh00 && saturation[i][j]<=0.05) {
	  //~ if(grayi[i][j]>=thresh00) {
	  if( gray[i][j]>=thresh00) {
	    white[i][j] = true;
	  }
	}

      // fill flood to get candidate
      for(int i=0; i<H; ++i)
	for(int j=0; j<W; ++j)
	  if(white[i][j] && !visited[i][j]) {
	    Frame f(i,i,j,j);
	    int counter = 0;
	    browse(i, j, visited, white, f, counter);
	    int h = f.get_height();
	    int w = f.get_width();
	    bool valid = true;
	    //~ if(counter>20)
	      //~ cerr << "WHITESPOT C : " << f.y0 + h/2 << " " << f.x0 + w/2 << endl;
	    if(	valid && counter >=w*h/2 &&
		abs(h-w)<=5 && //3 
		MIN(h,w) >= 7 &&
		MAX(h,w) <= 15 ) {
	      int cy = f.y0 + h/2;
	      int cx = f.x0 + w/2;
	      // gradient intensity variance
	      double gradi, gradv;
	      {
		int shift=3;
		Frame fs( MAX(0, f.y0-shift), 
			  MIN(H, f.y1+1+shift), 
			  MAX(0, f.x0-shift), 
			  MIN(W, f.x1+1+shift) );
		int hs = fs.get_height();
		int ws = fs.get_width();
		gradi = cellsum(cumgrad,fs)/(double)(ws*hs);
		gradv = cellsum(cumgrad2,fs)/(double)(ws*hs);
		gradv -= gradi*gradi;
	      }	      
	      // Fixed threshold
	      double thresh1 = 0.7;
	      double thresh2 = 0.4;
	      if( gradi>thresh1 && gradv>thresh2 ) {
		cy = MAX(30, MIN(H-31, cy )); 
		cx = MAX(30, MIN(W-31, cx ));  
		Special spe(cy, cx, 0, h*w);
		spe.f.init(f);
		special.push_back(spe);
	      }
	  }
	}
    }
    
    void find_darkblob()
    {
      bool investigate = false;
      int histo0[256];
      memset( histo0, 0, sizeof(histo0) );
      for(int i=0; i<H; ++i)
      for(int j=0; j<W; ++j) {
	int k=MIN(255, (int)grayi[i][j] );
	++histo0[k];
      }
      int histo[256];
      memset( histo, 0, sizeof(histo) );
      // scharr like smoothing
      for(int i=1; i<255; ++i)
	histo[i] = 2*histo0[i] + histo0[i-1] + histo0[i+1];
      
      int* pmax = max_element(histo, histo+256);
      int imax = (int)(pmax - histo);
      int lmin = *pmax;
      int imin = imax;
      int imax2 = 0;
      for(int i=imax; i>=0; --i) {
	if(histo[i]<lmin)
	  lmin = histo[i], imin = i;
	if(histo[i]>lmin*2 && (lmin>100 || histo[i]>200) ) {
	  imax2 = i;
	  investigate = true;
	  break;
	}
      }
      int threshold = (imax2 + imin)/2;
      
      if(investigate) { 
	
	vvb dark(H, vb(W, false));
	vvb visited(H, vb(W, false));      
	for(int i=0; i<H; ++i)
	  for(int j=0; j<W; ++j) {
	    if(grayi[i][j]<=threshold) {
	      dark[i][j] = true;
	    }
	  }
	
	// fill flood to get candidate
	for(int i=0; i<H; ++i)
	  for(int j=0; j<W; ++j)
	    if(dark[i][j] && !visited[i][j]) {
	      Frame f(i,i,j,j);
	      int counter = 0;
	      browse(i, j, visited, dark, f, counter);
	      if(counter>1000) {
		badb.push_back(BadBlob(f));
	      }
	    }
	}
    }

    void find_whiteblob()
    {
      bool investigate = false;
      int histo0[256];
      memset( histo0, 0, sizeof(histo0) );
      for(int i=0; i<H; ++i)
      for(int j=0; j<W; ++j) {
	int k=MIN(255, (int)grayi[i][j] );
	++histo0[k];
      }
      int histo[256];
      memset( histo, 0, sizeof(histo) );
      for(int i=1; i<255; ++i)
	histo[i] = 2*histo0[i] + histo0[i-1] + histo0[i+1];
      
      int* pmax = max_element(histo, histo+256);
      int imax = (int)(pmax - histo);
      int lmin = *pmax;
      for(int i=imax; i<256; ++i) {
	if(histo[i]<lmin)
	  lmin = histo[i];
	if(i>=250 && (histo[i]>200) ) {
	  investigate = true;
	  break;
	}
      }
      investigate = true;
      int threshold = 250;
      if(investigate) { 
	
	vvb dark(H, vb(W, false));
	vvb visited(H, vb(W, false));      
	for(int i=0; i<H; ++i)
	  for(int j=0; j<W; ++j) {
	    if(grayi[i][j]>=threshold) {
	      dark[i][j] = true;
	    }
	  }
	
	// fill flood to get candidate
	for(int i=0; i<H; ++i)
	  for(int j=0; j<W; ++j)
	    if(dark[i][j] && !visited[i][j]) {
	      Frame f(i,i,j,j);
	      int counter = 0;
	      browse(i, j, visited, dark, f, counter);
	      if(counter>1000) {
		badb.push_back(BadBlob(f));
	      }
	    }
	}
    }

    void build_hsv()
    {
      #ifdef VERBOSE
	cerr << "start build_hsv" << endl;
      #endif
      hue.assign(H, vd(W, 0.0));
      value.assign(H, vd(W, 0.0));
      saturation.assign(H, vd(W, 0.0));
      double factor = 1.0/255.0;
      for(int i=0; i<H; ++i) {
	for(int j=0; j<W; ++j) {
	  int rr = rdata[i][j];
	  int gg = gdata[i][j];
	  int bb = bdata[i][j];
	  double fr = (double)rr*factor;
	  double fg = (double)gg*factor;
	  double fb = (double)bb*factor;
	  double cmin = (fr<fg) ? ( (fr<fb) ? fr : fb ) :  ( (fg<fb) ? fg : fb );
	  double cmax = (fr>fg) ? ( (fr>fb) ? fr : fb ) :  ( (fg>fb) ? fg : fb );
	  double delta = cmax - cmin;
	  // hue
	  if(delta<1e-6)
	    hue[i][j] = 0.0;
	  else if(rr>=gg && rr>=bb) {
	    hue[i][j] = 60*fmodf((fg-fb)/delta,6.0);
	  } else if( gg>=bb) {
	    hue[i][j] = 60*( (fb-fr)/delta + 2.0 );
	  } else {
	    hue[i][j] = 60*( (fr-fg)/delta + 4.0 );
	  }
	  // saturation
	  if(delta<1e-6)
	    saturation[i][j] = 0.0;
	  else {
	    saturation[i][j] = delta / cmax;
	  }
	  // value
	  value[i][j] = cmax;
	}
      }
    }
    
    void build_prov() {
      #ifdef VERBOSE
	cerr << "start build_prov" << endl;
      #endif
      build_gray();
      build_grad();
      normalize_color();
      build_hsv();
    }
    
    void clean_prov() {
      #ifdef VERBOSE
	cerr << "start clean_prov" << endl;
      #endif
      clear_vv(rdata);
      clear_vv(gdata);
      clear_vv(bdata);
      clear_vv(gray);
      clear_vv(grayi);
      clear_vv(grad);
      clear_vv(hue);
      clear_vv(saturation);
      clear_vv(value);
      clear_vv(cumgrad);
      clear_vv(cumgrad2);
      clear_vv(harris);
    }
};


class StructureRecognition
{
  int annotations_N;
  vector<Image> img;
  Image live;
  mapi ix;
  double tend_retrieval;
  int training_pos;
  int vote_modern;
  int vote_other;
  Grid grid0;
  
  public:
  string init(vs imgdata, int n)
  {
    #ifdef VERBOSE
      cerr << "start init" << endl;
    #endif
    initTime();
    // limit of the training 
    tend_retrieval = 60.0*20.0;
    grid0.set(30, 30, 10, 10);
    modern = true;
    vote_modern = 0;
    vote_other = 0;
    // Max number of annotations allowed per testcase = 500 .. 8000;    
    annotations_N = n;
    img.clear();
    img.reserve(imgdata.size());
    
    {
      int i=0;
      for(auto s=imgdata.begin(), send=imgdata.end(); s<send; ++s,++i) {
	//parse
	int pos=0;
	int k= count(s->begin(), s->end(), ',');
	k /= 3;
	for(;(*s)[pos]!=',';++pos);
	Image newimg;
	newimg.id = s->substr(0, pos);
	ix[newimg.id] = i;
	
	for(int j=0; j<k; ++j) {
	  int nbefore = atoi( &(*s)[++pos] );
	  for(;(*s)[pos]!=',';++pos);
	  int nlifetime = atoi( &(*s)[++pos] );
	  for(;(*s)[pos]!=',';++pos);
	  int time = atoi( &(*s)[++pos] );
	  newimg.annotator.push_back( Annotator( nbefore, nlifetime, time) );
	  if(j!=k-1)
	    for(;(*s)[pos]!=',';++pos);
	}
	img.push_back(newimg);
      }
    }
    
    training_pos = 0;
    if(skip_training) {
      modern = true;
      return "END";
    } else {
      return img[training_pos].id;
    }
  }
  
  string end_of_training() 
  {
    modern = (vote_modern>2*vote_other); 
    cerr << "modern/other :" << vote_modern << " / " << vote_other << endl; 
    // no time left
    return "END";    
  }
    
  string receiveImage(string imgid, vi imgdata)
  {
    #ifdef VERBOSE
      cerr << "start receiveImage" << endl;
    #endif
    // get image
    int idx = ix[imgid];
    Image& cur = img[idx];
    cur.receive_data(imgdata);
    cur.build_prov();
    cur.sinit();
    cur.score_training(grid0);
    cur.clean_prov();
    
    // if not out of time, get 1st annotation 
    // decide next action
    if( cur.call_for_annotations )
      return ( cur.id + ",0" );
    else if(Time()>=tend_retrieval) { // no time left
      return end_of_training();
    } else if (training_pos < img.size()-1) {// get next image
      return img[++training_pos].id;
    }
    return end_of_training();
  }

  
  string receiveAnnotations(string imgid, int a, vs annotations)
  {
    #ifdef VERBOSE
      cerr << "start receiveAnnotations" << endl;
    #endif
    int idx = ix[imgid];
    Image& cur = img[idx];
    // parse annotation
    for(auto s=annotations.begin(), send=annotations.end(); s<send; ++s) {
      int pos=0;
      for(;(*s)[pos]!=',';++pos);
      string atype = s->substr(0, pos);
      int x = atoi( &(*s)[++pos] );
      for(;(*s)[pos]!=',';++pos);
      int y = atoi( &(*s)[++pos] );
      if(atype=="other")
	cur.annotator[a].other.push_back(Point(y,x));
      else
	cur.annotator[a].target.push_back(Point(y,x));
    }
    
    // Cannot call annotations anymore : end training
    if(annotations_N<8) {
      cur.update_score_spe();
      vote_modern += cur.vote_modern;
      vote_other += cur.vote_other;
      return end_of_training();
    }
    
    // Request next annotation or next image
    --annotations_N;
    
    ++a;
    int n = cur.annotator.size();
    // Arbitrary max number of annotation retrieved per annotator
    // to avoid exhausting annotations_N on few images.
    if( a<MIN(25, n) )
      return ( cur.id + "," + int2string(a) );
    // no more annotation so check results for this image
    cur.update_score_spe();
    vote_modern += cur.vote_modern;
    vote_other += cur.vote_other;
    cerr << "Annotations credit : " <<  annotations_N << ", at " << imgid << endl;
    cerr << "Modern/other :" << vote_modern << " / " << vote_other << endl;
    printRam();
    cerr << " , Time : " << Time() << endl;

    if( Time()>=tend_retrieval || MAX(vote_modern, vote_other)>150 ) {
      // no time left
      return end_of_training();
    } else if (training_pos < img.size()-1) {
      // get next image
      return img[++training_pos].id;
    }
    // no image left
    return end_of_training();
  }
  
  vs labelImage(vi imgdata)
  {
    vs ret;
    #ifdef VERBOSE
      cerr << "start labelImage" << endl;
    #endif
    // get image
    Image& cur = live;
    cur.receive_data(imgdata);
    cur.build_prov();
    cur.sinit();
    cur.score_live(grid0);
    cur.clean_prov();
    for(auto s=cur.label.begin(), send=cur.label.end(); s!=send; ++s) {
      stringstream ss;
      ss << s->x << "," << s->y << "," << s->proba ;
      ret.push_back(ss.str());
    }
    return ret;
  }
};

// main method for interaction with the visualizer
// added by TopCoder
int main() {
	int N, len;
	cin >> N >> len;

	vector<string> metaData(len);
	for (int i=0; i < len; i++) {
		cin >> metaData[i];
	}

	StructureRecognition obj;
	string cmd = obj.init(metaData, N);

	cout << cmd << endl << flush;

	while (1) {
		int visualcmd;
		cin >> visualcmd;

		if (visualcmd == 0) {
			break;
		}

		if (visualcmd == 1) {
			string imageId;
			int imageSize;
			cin >> imageId >> imageSize;
			vector<int> imageData(imageSize);
			for (int i=0; i < imageSize; i++) {
				cin >> imageData[i];
			}
			string cmd = obj.receiveImage(imageId, imageData);
			cout << cmd << endl << flush;
		}

		if (visualcmd == 2) {
			string imageId;
			int annotatorIndex, annotatorLen;
			cin >> imageId >> annotatorIndex >> annotatorLen;
			vector<string> annotations(annotatorLen);
			for (int i=0; i < annotatorLen; i++) {
				cin >> annotations[i];
			}
			string cmd = obj.receiveAnnotations(imageId, annotatorIndex, annotations);
			cout << cmd << endl << flush;
		}

		if (visualcmd == 3) {
			int imageSize;
			cin >> imageSize;
			vector<int> imageData(imageSize);
			for (int i=0; i < imageSize; i++) {
				cin >> imageData[i];
			}
			vector<string> labels = obj.labelImage(imageData);
			cout << labels.size() << endl;
			for (int i=0; i < labels.size(); i++) {
				cout << labels[i] << endl;
			}
			cout << flush;
		}
	}

	return 0;
}