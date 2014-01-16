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

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <ctype.h>
#include <time.h>
#include <assert.h>

#include <vector>
#include <list>
#include <map>
#include <set>
#include <queue>
#include <deque>
#include <stack>
#include <bitset>
#include <algorithm>
#include <functional>
#include <numeric>
#include <utility>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <complex>
#include <fstream>

using namespace std;

#define FOR(i, b, e)    for(int i = (b); i <= (e); i++)
#define FORL(i, b, e)    for(int i = (b); i < (e); i++)
#define FORD(i, e, b)    for(int i = (e); i >= (b); i--)
#define FOR0(i, e)        FORL(i, 0, e)

#define min(a, b)        (((a) < (b)) ? (a) : (b))
#define max(a, b)        (((a) > (b)) ? (a) : (b))
#define MINA(a, b)        do { if ((a) > (b)) (a) = (b); } while(0)
#define MAXA(a, b)        do { if ((a) < (b)) (a) = (b); } while(0)
#define MINA2(a, b, i, j)        do { if ((a) > (b)) { (a) = (b); (i) = (j); } } while(0)
#define MAXA2(a, b, i, j)        do { if ((a) < (b)) { (a) = (b); (i) = (j); } } while(0)

#define SWAP(a, b)        do { int _t = a; a = b; b = _t; } while(0)
#define SWAPT(a, b, t)    do { t _t = a; a = b; b = _t; } while(0)
#define SQR(a)            ((a) * (a))
#define MSET(a, b)        memset(a, b, sizeof(a))

#define INT                int
#define INT_CAP            0x3F3F3F3F

typedef long long int    LI;

typedef pair<int, int>    II;
typedef vector<int>       VI;
typedef vector<double>    VD;
typedef vector<string>    VS;
#define ALL(c)            c.begin(), c.end()
#define SZ(c)            (static_cast<int>(c.size()))
#define FORALL(it, c)     for(it = c.begin(); it != c.end(); ++it)
#define FORALLR(it, c)    for(it = c.rbegin(); it != c.rend(); ++it)
#define PB                push_back
#define MP                make_pair
#define P1                first
#define P2                second

FILE *debugf = stdout;
char *testcase = "";
char *testmode = "";
int tester, testdump, testparam, testcommand;
int timerstep, timerstepcap, localtimeout;
LI globalstart, localstart;

#define DEBUG(f, ...)    do { fprintf(debugf, "DEBUG:" f "\n", ##__VA_ARGS__); fflush(debugf); } while(0) 

#ifdef _DEBUG
#define DEBUG1(f, ...)    do { fprintf(debugf, "DEBUG:" f "\n", ##__VA_ARGS__); fflush(debugf); } while(0) 
#define ASSERT	assert
#else
#define DEBUG1(f, ...)	  do {} while(0)
#define ASSERT
#endif

#ifndef LOCAL

#define DEBUGOUT	cout
#include <sys/time.h>
 
LI gettimems() 
{
    timeval tv; 
    gettimeofday(&tv, 0); 
    return tv.tv_sec * 1000 + tv.tv_usec / 1000;
}
#else

#define DEBUGOUT	cerr
#define sprintf sprintf_s

#include <time.h>

LI gettimems() 
{
    return (LI) clock() * 1000 / CLOCKS_PER_SEC;
}

#endif

int randn0(int n)
{
/*	DEBUG("RANDDDD!");
	if (n <= 1)
		return 0; */
//	return 1;  
	return rand() % n; 				/* fast, but incorrect */
}

string itos(int i)
{
	ostringstream os;
	os << i;
	return os.str();
}

string dbltos(double dbl)
{
	ostringstream os;
	os << dbl;
	return os.str();
}

unsigned int i32sqrt(unsigned int n)   
{   
	unsigned int c, g;
	c = g = 0x8000;   
	while(1) {   
		if ((unsigned long long) g * g > n)   
			g ^= c;   
		c >>= 1;   
		if (!c)   
			return g;   
		g |= c;   
	}   
}  

/*---*/

#define COORD_DOUBLE
#define GETDBL(x)		x = COORD_CAP

#define PI		3.141592653589793L
#define RAD(x)	(PI * (x) / 180)
#define EPS		1e-4 // must be set
#define ROUND(x) ((x)>=0?(int)((x)+0.5):(int)((x)-0.5))
#define DBL_CAP	(DBL_MAX / 2)
#define INT_CAP	0x3F3F3F3F

//----------------------------------

#define MAX_W	1280
#define MAX_H	768

#define MINIHW	80
#define MINIHW2	(MINIHW + 2)

#define MAX_IMG			1000
#define MAX_IMG_ANNOTBLOCK	256
#define MAX_REL_ANNOT	5
#define MAX_IMG_ANNOT	(MAX_IMG_ANNOTBLOCK * MAX_REL_ANNOT)
#define MAX_INTENSITY	768
#define MAX_COLORVAL	8
#define MAX_HIT			256
#define MAX_PATTERN		6

#ifdef LOCAL
#define TIMELIMIT		35500000
#define TIMELIMIT1		35500000
#define TRAIN_TIMELIMIT	18000000
#else
#define TIMELIMIT		3540000
#define TIMELIMIT1		3040000
#define TRAIN_TIMELIMIT	1200000
#endif

enum {
	MODE_SINGLE = 0,
	MODE_CHECKMINI = 1,
	MODE_HITSTAT = 2,
	MODE_REAL = 3
};

enum {
	TARGET_MIX = 0,
	TARGET_MODERN = 1,
	TARGET_ANCIENT = 2
};

enum {
	HIT_NONE = -1,
	HIT_RIVER = 0,
	HIT_WHITETENT = 1,
	HIT_DARKRUIN = 2,
	HIT_LIGHTRUIN = 3,
	HIT_BWRUIN = 4,
	HIT_ROAD = 5,
	HIT_BLACKEDGE = 6,
	HIT_REDWHITEROOF = 7,
	HIT_BLUEROOF = 8,
	HIT_DARKSHAPE = 9,
	HIT_WHITERUIN = 10,
	HIT_MIXRUIN = 11,
	HIT_GREYRUIN = 12,
	HIT_BLACKRUIN = 13,
	HIT_DARKROOF = 14,
	HIT_SMALLRUIN = 15,
	HIT_INVERSERUIN = 16,
	HIT_3RUIN = 17,
	MAX_HITTYPE
};

string hitstrs[] = {"R", "wtent", "dr", "lr", "bwr", "r", "be", "rwroof", "broof", "ds", "wr", "mr", "gr", "br", "droof", "sr", "ir", "3r"};

enum {
	CATEGORY_MODERN = 0,
	CATEGORY_RIVER,
	CATEGORY_ROAD,
	CATEGORY_ANCIENT,
	CATEGORY_NONE,
	MAX_CATEGORY
};

char categories[MAX_CATEGORY] = {'M', 'R', 'r', 'a', 'n'};

struct annotblock_t {
	int curr;
	int total;
	int sec;
};

struct annot_t {
	int annotblock;
	int x;
	int y;
	int sumx;
	int sumy;
	int structure;
	int other;
};

struct hit_t {
	int x;
	int y;
	int w;
	int h;
	int type;
	int cnt;
	int bestcnt;
	double err;
	double prob;
	int keyx;
	int mergecnt;
	int nroadriver;
};

struct pos_t {
	int x;
	int y;
};

struct catstat_t {
	int good;
	int badmodern;
	int badancient;
	int roadriver;
	int none;
};

struct hitstat_t {
	double goodrate;
	double scorerate;
	int good;
	int bad;
	int goodannot;
	int totalannot;
};

struct pattern_t {
	bool *quanttypes;
	int type;
};

struct annotspec_t {
	char category;
	int x;
	int y;
};

pos_t dirs[8] = {{1, 0}, {1, 1}, {0, 1}, {-1, 1}, {-1, 0}, {-1, -1}, {-1, 0}, {-1, 1}};

double quantrates[MAX_COLORVAL] = {0.005, 0.1, 0.785, 0.1, 1};
int quantvalues[MAX_COLORVAL] = {0, 0x404040, 0x808080, 0xC0C0C0, 0xFFFFFF, 0xFF0000, 0x0000FF};

int miniimg[MINIHW2 * MINIHW2];
int flood[MINIHW2 * MINIHW2];
int minicnt[(MINIHW + 1) * (MINIHW + 1)][8];

int img[MAX_H * MAX_W];
int saveimg[MAX_H * MAX_W];

int w, h, maxannotreq, nannotblock, annotblock, hw, hwsize, test, ntrain, train, ndistannot, annotreq;
int nroadriver, nancient, nmodern, nhit, nproc;
bool timewarn, traintimewarn, timewarn1;
LI accumulatedtime;
int mode = MODE_REAL;
int target = TARGET_MIX;

string name;
int nannotblocks[MAX_IMG];
string ids[MAX_IMG];
hit_t hits[MAX_HIT];

annotblock_t annotblocks[MAX_IMG][MAX_IMG_ANNOTBLOCK];
annot_t distannots[MAX_IMG_ANNOT];
catstat_t catstats[MAX_CATEGORY];

hitstat_t hitstats[MAX_HITTYPE] = {{0}, {0.79,0.5}, {0.14,0.25}, {0.07,0.19}, {0},
{0}, {0.24, 0.32}, {0.3, 0.3}, {0.4, 0.5}, {0.4, 0.27},
{0.1, 0.2}, {0.15, 0.25}, {0.2, 0.25}, {0.15, 0.3}, {0.06, 0.28},
{0.11, 0.36}, {0.1, 0.2}, {0.1, 0.2}};
// precalculated ancient local rates are based on multiple hits, so actual goodrate is better, scorerate is weaker
// also, as it include hits from "other" category as well, scorerate is weakened with * 0.75

annotspec_t annotspecs[MAX_IMG][MAX_IMG_ANNOTBLOCK];
int nannotspecs[MAX_IMG];

bool blackonly[MAX_COLORVAL] = {1};
bool whiteonly[MAX_COLORVAL] = {0, 0, 0, 0, 1};
bool dark[MAX_COLORVAL] = {1, 1};
bool darkgreyonly[MAX_COLORVAL] = {0, 1};
bool midgreyonly[MAX_COLORVAL] = {0, 0, 1};
bool light[MAX_COLORVAL] = {0, 0, 0, 1, 1};
bool blackwhite[MAX_COLORVAL] = {1, 0, 0, 0, 1};
bool redwhite[MAX_COLORVAL] = {0, 0, 0, 0, 1, 1};
bool blueonly[MAX_COLORVAL] = {0, 0, 0, 0, 0, 0, 1};
bool mix[MAX_COLORVAL] = {1, 1, 0, 1, 1, 1, 1};

#ifdef LOCAL
ofstream newannotf;
#endif

pattern_t patterns[MAX_PATTERN] = {{blackonly, HIT_BLACKRUIN}, {whiteonly, HIT_WHITERUIN}, {dark, HIT_DARKRUIN}, {light, HIT_LIGHTRUIN}, {darkgreyonly, HIT_GREYRUIN}, {mix, HIT_MIXRUIN}};

int circlerowws[MINIHW / 2 + 1][MINIHW];

int* minip(int r, int c)
{
	return &miniimg[(r + 1) * MINIHW2 + c + 1];
}

int getminip(int r, int c)
{
	return *minip(r, c);
}

void setminip(int r, int c, int v)
{
	*minip(r, c) = v;
}

int* imgp(int r, int c)
{
	return &img[r * w + c];
}

int getp(int r, int c)
{
	return *imgp(r, c);
}

void setp(int r, int c, int v)
{
	*imgp(r, c) = v;
}

int red(int rgb)
{
	return (rgb >> 16) & 0xFF;
}

int green(int rgb)
{
	return (rgb >> 8) & 0xFF;
}

int blue(int rgb)
{
	return rgb & 0xFF;
}

void sumup()
{
	memset(minicnt, 0, MINIHW * 8 * sizeof(int));
	FOR0(q, MAX_COLORVAL) {
		int qv = quantvalues[q];
		FOR0(r, MINIHW) {
			minicnt[(r + 1) * (MINIHW + 1)][q] = 0;
			FOR0(c, MINIHW) {
				int p = (r + 1) * (MINIHW + 1) + c + 1;
				minicnt[p][q] = minicnt[p - 1][q] + minicnt[r * (MINIHW + 1) + c + 1][q] - minicnt[r * (MINIHW + 1) + c][q];
				if (*minip(r, c) == qv)
					minicnt[p][q]++;
			}
		}
	}
}

int countblock(int bx, int by, int bw, int bh, int q)
{
	return minicnt[(by + bh) * (MINIHW + 1) + bx + bw][q] - minicnt[(by + bh) * (MINIHW + 1) + bx][q] -
		minicnt[by * (MINIHW + 1) + bx + bw][q] + minicnt[by * (MINIHW + 1) + bx][q];
}

int countblock(int bx, int by, int bw, int bh, bool* quanttypes)
{
	int cnt = 0;
	FOR0(q, MAX_COLORVAL) {
		if (!quanttypes[q])
			continue;
			cnt += countblock(bx, by, bw, bh, q);
	}
	return cnt;
}

int countblockaround(int bx, int by, int bw, int bh, bool* quanttypes, int d, int& totalcnt)
{
	int x = max(0, bx - d);
	int y = max(0, by - d);
	int w = min(MINIHW - x, bw + 2 * d);
	int h = min(MINIHW - y, bh + 2 * d);
	totalcnt = h * w;
	int cnt = countblock(x, y, w, h, quanttypes);
	return cnt;
}

int countcircle(int bx, int by, int r, int q, int halfmode)
{
	int cnt = 0;
	int rend = 2 * r;
	if (halfmode)
		rend = r;
	FOR0(dy, rend) {
		int w = halfmode == 2 ? circlerowws[r][r - dy] : circlerowws[r][dy];
		cnt += countblock(bx + r - w / 2, by + dy, w, 1, q);
	}
	return cnt;
}

int countcircle(int bx, int by, int r, bool quanttypes[MAX_COLORVAL], int halfmode)
{
	int cnt = 0;
	FOR0(q, MAX_COLORVAL) {
		if (!quanttypes[q])
			continue;
		cnt += countcircle(bx, by, r, q, halfmode);
	}
	return cnt;
}

bool checkblock(int bx, int by, int bw, int bh, bool quanttypes[MAX_COLORVAL], double minrate, double maxrate, double& rate)
{
	int bwh = bw * bh;
	int cnt = countblock(bx, by, bw, bh, quanttypes);
	rate = (double) cnt / bwh;
	return minrate <= rate && rate <= maxrate;
}

bool checkring(int bx, int by, int r, bool quanttypes[MAX_COLORVAL], double maxerr, double& err, int& incircle, int halfmode)
{
	int bwh = SQR(2 * r);
	int bwcbig = (int) (SQR(r + 1) * PI);
	int bwcsmall = (int) (SQR(r - 1) * PI);
	int bwcring = bwcbig - bwcsmall;
	int inblock = countblockaround(bx, by, 2 * r, halfmode ? r : 2 * r, quanttypes, 1, bwh);
	if (inblock < (double) bwcring * (1.0 - maxerr)) { // testing that there are not many outside is unimportant
		err = maxerr + EPS;
		return false;
	}
	int inbig = countcircle(bx, by, r + 1, quanttypes, halfmode);
	int insmall = countcircle(bx, by, r - 2, quanttypes, halfmode);
	int inring = inbig - insmall;
	err = (double) (bwcring - inring) / bwcring;
	return err <= maxerr;
}


bool checkcircle(int bx, int by, int r, bool quanttypes[MAX_COLORVAL], double maxerr, double& err, int& incircle, int halfmode)
{
	int bwh = SQR(2 * r);
	int bwc = (int) (SQR(r) * PI);
	if (halfmode != 0) {
		bwh /= 2;
		bwc /= 2;
	}
	int inblock = countblockaround(bx, by, 2 * r, halfmode ? r : 2 * r, quanttypes, 1, bwh);
	int bwo = bwh - bwc;
	if (inblock < (double) bwc * (1.0 - maxerr) || 
		(inblock > bwc && inblock - bwc > (double) bwo * maxerr)) { // testing that there are not many outside is unimportant
		err = maxerr + EPS;
		return false;
	}
	incircle = countcircle(bx, by, r, quanttypes, halfmode);
	err = (double) (bwc - incircle) / bwc;
	err += (double) (inblock - incircle) / bwo;
	if (halfmode)
		err *= 1.5;
	return err <= maxerr;
}

void initcircleroww()
{
	FOR(r, 1, MINIHW / 2) {
		int r2 = SQR(r);
		FOR0(yd, 2 * r)
			circlerowws[r][yd] = 2 * (int) i32sqrt(r2 - SQR(r - yd));
	}
}

#ifdef LOCAL

void savetga(string sfn, int w, int h, int lw, int *img)
{
	unsigned char buf[4 * MAX_W];
	ofstream ofs(sfn.c_str(), ios::out | ios::binary);
	ofs.put(0);
	ofs.put(0);
	ofs.put(2);
	FOR0(i, 9)
		ofs.put(0); 		
	ofs.put(w & 0xFF);
	ofs.put(w / 256);
	ofs.put(h & 0xFF);
	ofs.put(h / 256);
	ofs.put(32);
	ofs.put(0);
	FORD(y, h -1, 0) {
		unsigned char *p = buf;
		FOR0(x, w) {
			int rgb = img[y * lw + x];
			*p++ = blue(rgb);
			*p++ = green(rgb);
			*p++ = red(rgb);
			*p++ = 0;
		}
		ofs.write((char *) buf, 4 * w);
	}
}

void savetga(string sfn)
{
	savetga(sfn, w, h, w, img);
}

void savedat(string sfn)
{
	ofstream ofs(sfn.c_str(), ios::binary);
	ofs.write((char *)&h, sizeof(h));
	ofs.write((char *)&w, sizeof(w));
	ofs.write((char *)img, hwsize);
}

void loaddat(string sfn)
{
	ifstream ifs(sfn.c_str(), ios::binary);
	ifs.read((char *)&h, sizeof(h));
	ifs.read((char *)&w, sizeof(w));
	hw = h * w;
	hwsize = hw * sizeof(int);
	ifs.read((char *)img, hwsize);
}

void savemini(string sfn)
{
	savetga(sfn, MINIHW, MINIHW, MINIHW2, miniimg + MINIHW2 + 1);
}

#else

void savetga(string sfn, int w, int h, int lw, int *img)
{
}

void savetga(string sfn)
{
}

void savedat(string sfn)
{
}

void loaddat(string sfn)
{
}

void savemini(string sfn)
{
}

#endif

void pickmini(int x, int y)
{
	x = max(MINIHW / 2, min(w - MINIHW / 2, x));
	y = max(MINIHW / 2, min(h - MINIHW / 2, y));
	FOR0(r, MINIHW)
		memcpy(minip(r, 0), imgp(y - MINIHW / 2 + r, x - MINIHW / 2), MINIHW * sizeof(int));
}

bool checkcircles(bool quanttypes[MAX_COLORVAL], hit_t& best, double maxerr, int minr = 10, int maxr = 40, bool allowhalf = false, bool allowring = false)
{
	best.err = 1;
	best.bestcnt = 0;
	FORD(r, maxr, minr) {
		int halfmodecap = r >= 10 && allowhalf ? 2 : 0;
		int xstep = r > 15 ? 2 : 1;

		FOR(halfmode, 0, halfmodecap) {
			int ysize = halfmode ? r : 2 * r;
			int ystart = halfmode == 2 ? MINIHW - r - r / 4 : 0;
			int yend = halfmode == 1 ? r / 4 : MINIHW - ysize - 1;
			for(int y = ystart; y <= yend; y += 1) {
				int inblock = countblock(0, y, MINIHW, MINIHW - y, quanttypes);
				int bwh = SQR(2 * r);
				int bwc = (int) (SQR(r) * PI);
				int bwo = bwh - bwc;
				if (inblock < (double) bwc * (1.0 - maxerr)) // testing that there are not many outside is unimportant
					continue;
				for(int x = 0; x + 2 * r < MINIHW; x += xstep) {
					double err;
					int cnt;
//					if (r == 16 && x == 30 && y == 62 && halfmode == 1)
//						cnt = 0;
					if (checkcircle(x, y, r, quanttypes, maxerr, err, cnt, halfmode)) {
		//					if (err < 0.5)
		//						cerr << r << ' ' << x << ' ' << y << ' ' << err << endl;
						best.bestcnt++;
						if (err < best.err) {
							best.cnt = cnt;
							best.err = err;
							best.x = x;
							best.y = y;
							best.w = 2 * r;
							best.h = ysize;
						}
					}
				}
			}
		}
	}
	return best.err < 1;
}

bool checksetflood(int x, int y, int qv1, int qv2, bool& goodpix, int floodset, int floodunset)
{
	goodpix = false;
	int v = miniimg[y * MINIHW2 + x];
	if (v != qv1 && v != qv2)
		return false;
	goodpix = true;
	int f = flood[y * MINIHW2 + x];
	if (f != floodunset)
		return false;
	flood[y * MINIHW2 + x] = floodset;
	return true;
}

queue<pos_t> qq;

int doflood(int x, int y, int qv1, int qv2, bool diagonal, int& xmin, int& xmax, int& ymax, int& massivecnt, int floodset = 1, int floodunset = 0)
{
	x++;
	y++;
	bool goodpix;
	massivecnt = 0;
	if (!checksetflood(x, y, qv1, qv2, goodpix, floodset, floodunset))
		return 0;
	int cnt = 1;
	xmin = x;
	xmax = x;
	ymax = y;
	pos_t p;
	p.x = x;
	p.y = y;
	qq.push(p);
	while(!qq.empty()) {
		p = qq.front();
		qq.pop();
		int incnt = 0;
		for(int d = 0; d < 8; d += 1 + !diagonal) {
			pos_t p1;
			p1.x = p.x + dirs[d].x;
			p1.y = p.y + dirs[d].y;
			bool ok = checksetflood(p1.x, p1.y, qv1, qv2, goodpix, floodset, floodunset);
			if (goodpix)
				incnt++;
			if (!ok)
				continue;
			MINA(xmin, p1.x);
			MAXA(xmax, p1.x);
			MAXA(ymax, p1.y);
			qq.push(p1);
			cnt++;
		}
		if (incnt >= 4 + 3 * diagonal)
			massivecnt++;
	} 
	xmax--;
	xmin--;
	ymax--;
	return cnt;
}

enum {
	THICK_ANY = 0,
	THICK_MASSIVE = 1,
	THICK_THIN = 2,
	THICK_SUPERMASSIVE = 3
};

void clearflood()
{
	MSET(flood, 0);
	memset(flood, 0xFF, MINIHW2 * sizeof(int));
	memset(flood + (MINIHW + 1) * MINIHW2, 0xFF, MINIHW2 * sizeof(int));
	FOR(r, 1, MINIHW) {
		flood[r * MINIHW2] = -1;
		flood[r * MINIHW2 + MINIHW + 1] = -1;
	}
}

bool checkflood(bool quanttypes[MAX_COLORVAL], hit_t& best, int minhw = 5, int maxhw = MINIHW, int mode = THICK_MASSIVE, bool diagonal = false)
{
	best.cnt = 0;
	best.err = 0;
	best.bestcnt = 0;
	clearflood();
	int cap = minhw;
	if (mode == THICK_MASSIVE)
		cap = 4 * minhw;
	else if (mode == THICK_SUPERMASSIVE)
		cap = 6 * minhw;
	if (countblock(0, 0, MINIHW, MINIHW, quanttypes) < cap)
		return false;
	int qv1 = 0x11111111;
	int qv2 = 0x11111111; // speedup only two active quant
	FOR0(q, MAX_COLORVAL) {
		if (quanttypes[q]) {
			if (qv1 == qv2)
				qv1 = quantvalues[q];
			else
				qv2 = quantvalues[q];
		}
	}
	for(int y = 0; y < MINIHW; ++y) {
		for(int x = 0; x < MINIHW; ++x) {
			int xmin, xmax, ymax;
			int massivecnt;
			int cnt = doflood(x, y, qv1, qv2, diagonal, xmin, xmax, ymax, massivecnt);
			if (cnt > 10) {
				int w = xmax - xmin + 1;
				int h = ymax - y + 1;
				int mhw = max(h, w);
	//			if (searchmassive && cnt > 100)
	//				cnt = cnt;
				bool ok = true;
				switch(mode) {
				case THICK_MASSIVE:
					ok = massivecnt >= cnt / 2;
					break;
				case THICK_THIN:
					ok = 3 * massivecnt < cnt;
					break;
				case THICK_SUPERMASSIVE:
					ok = massivecnt > 3 * cnt / 4;
					break;
				}
				if (minhw <= mhw && mhw <= maxhw && ok) {
					best.bestcnt++;
					if (cnt > best.cnt) {
						best.cnt = cnt;
						best.x = xmin;
						best.y = y;
						best.w = w;
						best.h = h;
						best.keyx = x;
					}
				}
			}
		}
	}
	return best.cnt > 0;
}

bool rechecksoloblob(bool quanttypes[MAX_COLORVAL], hit_t& best, double maxerr)
{
	int rcnt = countblock(best.x, best.y, best.w, best.h, quanttypes);
	best.err = (double) rcnt / best.cnt - 1;
	return best.err <= maxerr;
}

bool hassomeblack(hit_t& best)
{
	int dummy = 0;
	int cnt = countblockaround(best.x, best.y, best.w, best.h, blackonly, 1, dummy);
	return cnt > (best.w + best.h) / 4 && cnt <= best.w + best.h;
}

int countlinearviews(bool quanttypes[MAX_COLORVAL], hit_t& best, int& maxdiff, bool diagonalfill = false, bool multiline = false)
{
	int qv1 = 0x11111111;
	int qv2 = 0x11111111; // speedup only two active quant
	FOR0(q, MAX_COLORVAL) {
		if (quanttypes[q]) {
			if (qv1 == qv2)
				qv1 = quantvalues[q];
			else
				qv2 = quantvalues[q];
		}
	}

	int dummy;
	doflood(best.keyx, best.y, qv1, qv2, diagonalfill, dummy, dummy, dummy, dummy, 2, 1);

	int firstdir, lastdir, prevfirst, prevlast, prevfirstdir, prevlastdir;
	int firstdir2, lastdir2, prevfirstdir2, prevlastdir2;
	int lastbad = 0;
	int firstbad = 0;
	maxdiff = 0;
	int cap = 0;
	int cnt = 0;
	int zerofirst = 0;
	int zerolast = 0;
	FOR(y, best.y, best.y + best.h - 1) {
		int first = -1;
		int last = -1;
		FOR(x, best.x, best.x + best.w - 1) {
			int v = flood[(y + 1) * MINIHW2 + x + 1];
			if (v == 2) {
				last = x;
				if (first == -1)
					first = x;
			}
		}
		MAXA(maxdiff, last - first);
		if (first == 0)
			zerofirst++;
		if (last == MINIHW - 1)
			zerolast++;
		if (y > best.y) {
			firstdir = first - prevfirst;
			lastdir = last - prevlast;
			if (y > best.y + 1) {
				firstdir2 = firstdir + prevfirstdir;
				lastdir2 = lastdir + prevlastdir;
				if (y > best.y + 2) {
					if (abs(firstdir2 - prevfirstdir2) > 1)
						firstbad++;
					if (abs(lastdir2 - prevlastdir2) > 1)
						lastbad++;
				}
				prevfirstdir2 = firstdir2;
				prevlastdir2 = lastdir2;
			}
			prevfirstdir = firstdir;
			prevlastdir = lastdir;
		}
		prevfirst = first;
		prevlast = last;
	}
	int badcap = best.h <= 10 ? 0 : 2 + best.h / 20;
	if ((multiline || maxdiff >= 10)) {
		if ((best.x + best.w < MINIHW || zerolast < best.h / 4) && lastbad <= badcap)
			cnt++;
		if ((best.x > 0 || zerofirst < best.h / 4) && firstbad <= badcap)
			cnt++;
	}

	lastbad = 0;
	firstbad = 0;
	zerofirst = 0;
	zerolast = 0;
	cap = best.w;
	FOR(x, best.x, best.x + best.w - 1) {
		int first = -1;
		int last = -1;
		FOR(y, best.y, best.y + best.h - 1) {
			int v = flood[(y + 1) * MINIHW2 + x + 1];
			if (v == 2) {
				last = y;
				if (first == -1)
					first = y;
			}
		}
		MAXA(maxdiff, last - first);
		if (first == 0)
			zerofirst++;
		if (last == MINIHW - 1)
			zerolast++;
		if (x > best.x) {
			firstdir = first - prevfirst;
			lastdir = last - prevlast;
			if (x > best.x + 1) {
				firstdir2 = firstdir + prevfirstdir;
				lastdir2 = lastdir + prevlastdir;
				if (x > best.x + 2) {
					if (abs(firstdir2 - prevfirstdir2) > 1)
						firstbad++;
					if (abs(lastdir2 - prevlastdir2) > 1)
						lastbad++;
				}
				prevfirstdir2 = firstdir2;
				prevlastdir2 = lastdir2;
			}
			prevfirstdir = firstdir;
			prevlastdir = lastdir;
		}
		prevfirst = first;
		prevlast = last;
	}
	badcap = best.w < 10 ? 0 : 2 + best.w / 20;
	if ((multiline || maxdiff >= 10)) {
		if ((best.y + best.h < MINIHW || zerolast < best.w / 4) && lastbad <= badcap)
			cnt++;
		if ((best.y > 0 || zerofirst < best.w / 4) && firstbad <= badcap)
			cnt++;
	}
	return cnt;
}

bool checkedges(bool quanttypes[MAX_COLORVAL], hit_t& best, bool diagonal = false)
{
	if (!checkflood(quanttypes, best, 20, MINIHW, THICK_THIN, diagonal))
		return false;

	if (best.bestcnt <= 1 && best.cnt <= 5 * max(best.h, best.w) / 2) // one edge-combo is not enough unless it is big
		return false;

	bool multiline = best.bestcnt >= 2;

	int maxdiff;
	int cntlin = countlinearviews(quanttypes, best, maxdiff, diagonal, multiline);
	
	if (cntlin < 2)
		return false;

	if (multiline && maxdiff < 10 && checkflood(light, best, MINIHW / 3, MINIHW, THICK_THIN, true))
		return false;

	return true;
}

bool checkroof(bool quanttypes[MAX_COLORVAL], hit_t& best, int minhw = 10)
{
	if (!checkflood(quanttypes, best, minhw))
		return false;
	int dummy;
	int cnt = countlinearviews(quanttypes, best, dummy);
	if (quanttypes[1] && cnt < 3) // grey should be more selective
		return false;
	if (cnt < 2)
		return false;
	if (quanttypes[4]) {
		if (!hassomeblack(best))
			return false;
		return rechecksoloblob(light, best, 0.3); // for white, we dont want lightgrey either
	} 
	return rechecksoloblob(quanttypes, best, 0.2);
}


bool checkall(hit_t& best)
{
//	int dummy;
	sumup();
	if (checkroof(blueonly, best)) {
		best.type = HIT_BLUEROOF;
		return true;
	}
	if (checkedges(blackonly, best)) {
		best.type = HIT_BLACKEDGE;
		return true;
	} 
	if (checkcircles(whiteonly, best, 0.46, 6, 10, false)) {
		hit_t dummy;
		if (best.err < 0.3 || hassomeblack(best) && !checkflood(whiteonly, dummy, MINIHW / 2)) { // either good hit or has some blacks as well
			best.type = HIT_WHITETENT;
			return true;
		}
	}

	if (checkflood(light, best, MINIHW)) {
		best.type = HIT_ROAD;
		return true;
	}
	if (checkflood(dark, best,  MINIHW)) {
		best.type = HIT_RIVER;
		return true;
	}
	if (checkroof(redwhite, best)) {
		best.type = HIT_REDWHITEROOF;
		return true;
	}

	if (target != TARGET_ANCIENT) {
		if (checkroof(darkgreyonly, best, 15)) {
			best.type = HIT_DARKROOF;
			return true;
		}
	}

	if (checkedges(blackonly, best, true)) { 
		if (target != TARGET_ANCIENT) { // skip it, as overlaps with dark ruin
			best.type = HIT_BLACKEDGE;
			best.err = 0.5;
			return true;
		} 
		else {
			best.type = HIT_DARKRUIN;
			return true;
		} 
	}

	if (target == TARGET_MODERN) {
		if (checkcircles(dark, best, 0.4, 15) && hassomeblack(best)) {
			best.type = HIT_DARKSHAPE;
			return true;
		}
	} 

	best.err = 1;
	FOR0(i, MAX_PATTERN) {
		hit_t hit;
		checkcircles(patterns[i].quanttypes, hit, target == TARGET_ANCIENT ? 0.52 : 0.48);
		if (best.err > hit.err) {
			best = hit;
			best.type = patterns[i].type;
		}
	}
	if (best.err < 1)
		return true;

	if (target == TARGET_ANCIENT) {
		if (checkcircles(darkgreyonly, best, 0.4, 5, 9) || 
			checkcircles(blackonly, best, 0.4, 5, 9) ||
			checkcircles(blackwhite, best, 0.35, 5, 9)) {
//			best.err += 0.12;
			best.type = HIT_SMALLRUIN;
			return true;
		} 
/*		if (checkcircles(darkgreyonly, best, 0.5, 2, 6) && best.bestcnt >= 3 && best.bestcnt <= 8) {
			best.type = HIT_3RUIN;
			return true;
		} */
	}

	if (checkflood(dark, best, 3 * MINIHW / 4, MINIHW, THICK_SUPERMASSIVE)) { // not so long, but very massive one
		best.type = HIT_RIVER;
		return true;
	}
	if (checkflood(light, best, 3 * MINIHW / 4, MINIHW, THICK_THIN)) { // thin light; may not be needed, as it doesnt make a difference
		best.type = HIT_ROAD; 
		return true;
	}

	if (checkflood(dark, best, 3 * MINIHW / 4, MINIHW, THICK_THIN)) {
		best.type = HIT_RIVER; 
		return true;
	}

	return false;
}

char resulttype(hit_t& best)
{
	switch (best.type) {
	case HIT_WHITETENT:
	case HIT_BLUEROOF:
	case HIT_BLACKEDGE:
	case HIT_REDWHITEROOF:
	case HIT_DARKSHAPE:
	case HIT_DARKROOF:
		return 'M';
	case HIT_RIVER:
		return 'R';
	case HIT_ROAD:
		return 'r';
	case HIT_NONE:
		return 'n';
	}
	return 'a';
}

void addtostat(hit_t& best, char expected)
{
	char result = resulttype(best);
	FOR0(i, MAX_CATEGORY) {
		if (categories[i] == expected) {
			if (tolower(expected) == tolower(result))
				catstats[i].good++;
			else if (result == 'M')
				catstats[i].badmodern++;
			else if (result == 'a')
				catstats[i].badancient++;
			else if (result == 'r' || result == 'R')
				catstats[i].roadriver++;
			else
				catstats[i].none++;
			break;
		}
	}
}

void writebest(hit_t& best)
{
#ifdef LOCAL
	char result = resulttype(best);
	cerr << result;
	if (result != 'n')
		cerr << ' ' << hitstrs[best.type] << ' ' << best.x << ' ' << best.y << ' ' << best.w << ' ' << best.h << ' ' << best.cnt << ' ' << best.err;
	cerr << endl; 
#endif
}

void calcbest(char expected)
{
	hit_t best;
	if (!checkall(best))
		best.type = HIT_NONE;
	addtostat(best, expected);
	writebest(best);
}

void init()
{
	initcircleroww();
}

bool contour(int x, int y, int dx, int dy)
{
	int v = getp(y, x);
	int sum = red(v) + green(v) + blue(v);
	int vn = getp(y - dy, x - dx);
	int sumn = red(vn) + green(vn) + blue(vn);
	if (sumn > 4 * sum / 3) {
		setp(y, x, (3 * red(v) / 4 << 16) + (3 * green(v) / 4 << 8) + 3 * blue(v) / 4); 
		return true;
	}
	return false;
}

void quant(int xd, int yd)
{
	int yr = h / yd;
	int xr = w / xd;
	int cap0 = hw / (xd * yd);
	FOR0(r0, yd) {
		FOR0(c0, xd) {
			int histo[MAX_INTENSITY];
			int lvl[MAX_INTENSITY];
			MSET(histo, 0);
			FOR(r, r0 * yr, (r0 + 1) * yr - 1) {
				FOR(c, c0 * xr, (c0 + 1) * xr - 1) {
					int v = getp(r, c);
					int sum = red(v) + green(v) + blue(v);
					histo[sum]++;
				}
			}
			int maxlvl = 0;
			int cnt = 0;
			int cap = (int) (cap0 * quantrates[maxlvl]);
			FOR0(i, MAX_INTENSITY) {
				if (cnt > cap) {
					maxlvl++;
					cap += (int) (cap0 * quantrates[maxlvl]);
				}
				lvl[i] = maxlvl;
				cnt += histo[i];
			}
			FOR(r, r0 * yr + 1, (r0 + 1) * yr - 1 - 1) {
				FOR(c, c0 * xr + 1, (c0 + 1) * xr - 1 - 1) {
					if (!contour(c, r, 1, 0))
						if (!contour(c, r, 0, 1))
							if (!contour(c, r, -1, 0))
								contour(c, r, 0, -1);
				}
			}
			FOR(r, r0 * yr, (r0 + 1) * yr - 1) {
				FOR(c, c0 * xr, (c0 + 1) * xr - 1) {
					int v = getp(r, c);
					int sum = red(v) + green(v) + blue(v);
					int nv;
					if (red(v) > 96 && 3 * red(v) >= 4 * green(v) && 3 * red(v) >= 4 * blue(v))
						nv = 0xFF0000;
					else if (blue(v) > 96 && 3 * blue(v) > 4 * red(v) && 3 * blue(v) >= 4 * green(v))
						nv = 0xFF;
					else
						nv = quantvalues[lvl[sum]];
					setp(r, c, nv);
				}
			}
		}
	}
}

void quantmini()
{
	int cap0 = SQR(MINIHW);
	int histo[MAX_INTENSITY];
	int lvl[MAX_INTENSITY];
	MSET(histo, 0);
	FOR0(r, MINIHW) {
		FOR0(c, MINIHW) {
			int v = getminip(r, c);
			int avg = (red(v) + green(v) + blue(v)) / 3;
			histo[avg]++;
		}
	}
	int maxlvl = 0;
	int cnt = 0;
	int cap = (int) (cap0 * quantrates[maxlvl]);
	FOR0(i, MAX_INTENSITY) {
		if (cnt > cap) {
			maxlvl++;
			cap += (int) (cap0 * quantrates[maxlvl]);
		}
		lvl[i] = maxlvl;
		cnt += histo[i];
	}
	FOR0(r, MINIHW) {
		FOR0(c, MINIHW) {
			int v = getminip(r, c);
			int avg = (red(v) + green(v) + blue(v)) / 3;
			int nv;
			nv = quantvalues[lvl[avg]];
			if (red(v) > 128 && 3 * red(v) > 4 * green(v) && 3 * red(v) > 4 * blue(v))
				nv = 0xFF0000;
			else if (blue(v) > 128 && blue(v) > 3 * (red(v) + green(v)) / 4)
				nv = 0xFF;
			setminip(r, c, nv);
		}
	}
}

void load(VI& image)
{
	h = image[0];
	w = image[1];
	hw = h * w;
	hwsize = hw * sizeof(int);
	int p = 2;
	FOR0(r, h)
		FOR0(c, w)
			setp(r, c, image[p++]);
}

void rebalance()
{
	int rsum = 0;
	int gsum = 0;
	int bsum = 0;
	FOR0(r, h) {
		FOR0(c, w) {
			int v = getp(r, c);
			rsum += red(v);
			gsum += green(v);
			bsum += blue(v);
		}
	}
	int tot = ((rsum + gsum + bsum) / 3) >> 8;
	rsum = ((rsum >> 8) + tot) / 2;
	gsum = ((gsum >> 8) + tot) / 2;
	bsum = ((bsum >> 8) + tot) / 2;
	FOR0(r, h) {
		FOR0(c, w) {
			int v = getp(r, c);
			v = (min(255, red(v) * tot / rsum) << 16) + (min(255, (green(v) * tot / gsum)) << 8) + min(255, blue(v) * tot / bsum);
			setp(r, c, v);
		}
	}
}

void preprocess()
{
	nproc = 0;
	savedat(name + ".dat");
	savetga(name + ".tga");
	rebalance();
	savetga(name + "_rb.tga");
	memcpy(saveimg, img, hwsize);
}

void doshrink()
{
	for(int r = 0; r < h; r += 2) {
		for(int c = 0; c < w; c += 2) {
			int rsum = 0;
			int gsum = 0;
			int bsum = 0;
			FOR0(rd, 2) {
				FOR0(cd, 2) {
					int v = getp(r + rd, c + cd);
					rsum += red(v);
					gsum += green(v);
					bsum += blue(v);
				}
			}
			FOR0(rd, 2) {
				FOR0(cd, 2) {
					setp(r + rd, c + cd, (rsum / 4 << 16) + (gsum / 4 << 8) + bsum / 4);
				}
			}
		}
	}
}

void process(int div, bool shrink = false)
{
	memcpy(img, saveimg, hwsize);
	if (shrink)
		doshrink();
	quant(div, div);
	string pref = 'q' + itos(div) + '_' + name;
	if (shrink)
		pref = 's'+ pref;
	savetga(pref + ".tga");
}

void processmeta(VS& meta)
{
	set<int> st;
	FOR0(i, SZ(meta)) {
		int pos = meta[i].find(',');
		ids[i] = meta[i].substr(0, pos);
//		cerr << "id: " << ids[i] << endl;
		while(pos != -1) {
			annotblock_t a;
			int pos1 = pos + 1;
			pos = meta[i].find(',', pos + 1);
			a.curr = atoi(meta[i].substr(pos1, pos - pos1).c_str());
			int pos2 = pos + 1;
			pos = meta[i].find(',', pos + 1);
			a.total = atoi(meta[i].substr(pos2, pos - pos2).c_str());
			pos2 = pos + 1;
			pos = meta[i].find(',', pos + 1);
			if (pos != -1)
				a.sec = atoi(meta[i].substr(pos2, pos - pos2).c_str());
			else
				a.sec = atoi(meta[i].substr(pos2).c_str());
			annotblocks[i][nannotblocks[i]] = a;
//			cerr << a.curr << ',' << a.total << ',' << a.sec << endl;
			nannotblocks[i]++;
		}
	}
//	cerr << "types:" << st.size() << endl;
//	for(set<int>::iterator it = st.begin(); it != st.end(); ++it)
//		cerr << *it << endl;
}

char get_annot_category(int x, int y)
{
	FOR0(i, nannotspecs[train]) {
		if (SQR(annotspecs[train][i].x - x) + SQR(annotspecs[train][i].y - y) <= SQR(40))
			return annotspecs[train][i].category;
	}
	return 0;
}

void addtoannot(int inserienr, int x, int y, int distannot)
{
#ifdef LOCAL
	newannotf << ',' << (get_annot_category(x, y) == 'a' ? '4' : '5') << ',' << x << ',' << y;
#endif
}

void processserie(VS& serie)
{
//	cerr << "serie on:" << SZ(serie) << endl;
#ifdef LOCAL
	newannotf << SZ(serie);
#endif
	FOR0(i, SZ(serie)) {
		int pos1 = serie[i].find(',') + 1;
		int pos = serie[i].find(',', pos1);
		int x = atoi(serie[i].substr(pos1, pos - pos1).c_str());
		int pos2 = pos + 1;
		pos = serie[i].find(',', pos + 1);
		int y = atoi(serie[i].substr(pos2, pos - pos2).c_str());
//		cerr << "serie:" << x << ',' << y << endl;
		bool structure = serie[i][0] == 's';
		bool repeat = false;
		if (mode == MODE_CHECKMINI || target == TARGET_ANCIENT) {
			FOR0(a, ndistannot) {
				if (SQR(distannots[a].x - x) + SQR(distannots[a].y - y) <= SQR(40)) { //! only merges closer hits with first pos, ok for now (may not be good for stat)
					addtoannot(i, x, y, a);
					distannots[a].sumx += x;
					distannots[a].sumy += y;
					if(structure)
						distannots[a].structure++;
					else
						distannots[a].other++;
					repeat = true;
					break;
				}
			}
			if (repeat)
				continue;
		}
		addtoannot(i, x, y, ndistannot);
		distannots[ndistannot].annotblock = annotblock;
		distannots[ndistannot].x = x;
		distannots[ndistannot].y = y;
		distannots[ndistannot].sumx = x;
		distannots[ndistannot].sumy = y;
		distannots[ndistannot].structure = structure;
		distannots[ndistannot].other = !structure;
		ndistannot++;
	}
#ifdef LOCAL
	newannotf << endl;
#endif
}

void finish()
{
	cerr << "finish" << endl;
#ifdef LOCAL
	newannotf.close();
#endif
	if (mode == MODE_HITSTAT || mode == MODE_REAL && target == TARGET_ANCIENT) {
#ifdef LOCAL
		ofstream ofs("mark.txt");
		ofs << train << endl;
		FOR0(i, MAX_HITTYPE)
			ofs <<  hitstats[i].good << ' ' << hitstats[i].bad << ' ' << hitstats[i].goodannot << ' ' << hitstats[i].totalannot << endl;
#endif
		FOR0(i, MAX_HITTYPE) {
			int cnt = hitstats[i].good + hitstats[i].bad;
			if (cnt > 10 && hitstats[i].good > 0 && hitstats[i].goodannot > 0) {
				hitstats[i].goodrate = (cnt ? (double) hitstats[i].good / cnt : 0);
				hitstats[i].scorerate = (hitstats[i].totalannot ? (double) hitstats[i].goodannot / hitstats[i].totalannot : 0);
#ifdef LOCAL
				if (target == TARGET_ANCIENT)
					hitstats[i].scorerate *= 0.75; // account for hits from other category
#endif
				cerr << "*";
			}
			cerr << hitstrs[i] << ' ' << hitstats[i].good << ' ' << hitstats[i].bad << ' ' << hitstats[i].goodrate;
			cerr << ' ' << hitstats[i].goodannot << ' ' << hitstats[i].totalannot << ' ' << hitstats[i].scorerate << endl;
		}
		return;
	}
	int allgood = 0;
	int allbadmodern = 0;
	int allbadancient = 0;
	int alltot = 0;
#ifdef LOCAL
	ofstream ofs("mark.txt");
	ofs << train << endl;
	FOR0(i, MAX_CATEGORY)
		ofs <<  catstats[i].good << ' ' << catstats[i].badmodern << ' ' << catstats[i].badancient << ' ' << catstats[i].roadriver << ' ' << catstats[i].none << endl;
#endif
	FOR0(i, MAX_CATEGORY) {
		int tot = catstats[i].good + catstats[i].badmodern + catstats[i].badancient + catstats[i].roadriver + catstats[i].none;
		alltot += tot;
		allbadmodern += catstats[i].badmodern;
		allbadancient += catstats[i].badancient;
		if (categories[i] == 'M' || categories[i] == 'a')
			allgood += catstats[i].good;
	}
	cerr <<  "cat" << '\t' << "good" << '\t' << "bad" << '\t' << "badmod" << '\t' << "badanc" << '\t' << "roadriv" << '\t' << "none" << endl;
	FOR0(i, MAX_CATEGORY) {
		int spec = categories[i] == 'M' ? allbadmodern : categories[i] == 'a' ? allbadancient : 0;
		cerr <<  categories[i] << '\t' << catstats[i].good << '\t' << spec << '\t' << catstats[i].badmodern << '\t' << catstats[i].badancient << '\t' << catstats[i].roadriver << '\t' << catstats[i].none << endl;
	}
	cerr << "total" << ' ' << alltot << ' ' << allgood << ' ' << allbadmodern << ' ' << allbadancient << ' ' << (double) allgood / alltot << ' ' << (double) (allbadmodern + allbadancient) / alltot << endl;
}

void processdistannots(bool check)
{
	cerr << "distannot:" << ndistannot << endl;
	FOR0(i, ndistannot) {
		if (mode == MODE_CHECKMINI && distannots[i].structure <= 2 && distannots[i].other <= 2)
			continue;
		if (mode == MODE_REAL && distannots[i].structure == 0) // for real run we don't care about nontarget
			continue;

		int cnt = distannots[i].structure + distannots[i].other;
		if (!check) {
			distannots[i].x = distannots[i].sumx / cnt;
			distannots[i].y = distannots[i].sumy / cnt;
		}

		string sfn;
		if (distannots[i].structure >= distannots[i].other)
			sfn = "modern_";
		else
			sfn = "other_";
		sfn += name + '_' + itos(distannots[i].annotblock) + '_' + itos(i);
		if (check)
			sfn += "_bw";
		sfn += ".tga";
		pickmini(distannots[i].x, distannots[i].y);
//		if (check)
//			quantmini();
		savemini(sfn);
		if (check) {
			char c = get_annot_category(distannots[i].x, distannots[i].y);
#ifdef LOCAL
			cerr << "data: " << train << ' ' << distannots[i].annotblock << ' ' << i << ' ' << distannots[i].x << ' ' << distannots[i].y << ' ' << distannots[i].structure << ' ' << distannots[i].other << ' ' << nannotblock << ' ' << c << ' ';
#endif
			if (!c) {
				if (distannots[i].structure >= distannots[i].other)
					c = 'M';
				else
					c = 'a';
			}
			calcbest(c);
		}
	}
	if (mode == MODE_REAL && target == TARGET_MIX) {
		cerr << "decision params " << catstats[CATEGORY_MODERN].good << ' ' << catstats[CATEGORY_MODERN].badancient << endl;
		if (catstats[CATEGORY_MODERN].good > 30 && catstats[CATEGORY_MODERN].good > 1.5 * catstats[CATEGORY_MODERN].badancient) {
			finish();
			cerr << "decision: modern" << endl;
			target = TARGET_MODERN;
			return;
		} else if (catstats[CATEGORY_MODERN].badancient > 60 && catstats[CATEGORY_MODERN].badancient > 3 * catstats[CATEGORY_ANCIENT].good) { // as we detect them as ancients
			finish();
			cerr << "decision: ancient" << endl;
			target = TARGET_ANCIENT;
			return;
		}
	}
//	ofstream ofs("mark.txt");
//	ofs << train;
}

void statdistannots()
{
	FOR0(h, nhit) {
		bool ok = false;
		int t = hits[h].type;
		FOR0(i, ndistannot) {
			int cnt = distannots[i].structure + distannots[i].other;
			int x = distannots[i].sumx / cnt;
			int y = distannots[i].sumy / cnt;
			char cat = 'n';
			if (distannots[i].structure) {
				if (target == TARGET_MODERN)
					cat = 'M';
				else
					cat = 'a';
			}
			bool goodtype = target == TARGET_MODERN && cat == 'M' ||
				target == TARGET_ANCIENT && cat == 'a';
			if (goodtype && SQR(x - hits[h].x) + SQR(y - hits[h].y) <= SQR(40)) {
				if (!(target == TARGET_MODERN && nhit > 6 || target == TARGET_ANCIENT && nhit > 10))
					hitstats[t].goodannot += distannots[i].structure;
				ok = true;
			}
		}
		if (!ok)
			hitstats[t].bad++;
		else {
			hitstats[t].good++;
			if (!(target == TARGET_MODERN && nhit > 6 || target == TARGET_ANCIENT && nhit > 10))
				hitstats[t].totalannot += nannotblock;
		}
	}
	finish();
}

void adjusttoedge(hit_t& best)
{
	if (best.x < 40)
		best.x = (best.x + 40) / 2;
	else if (best.x > 1196)
		best.x = (best.x + 1196) / 2;
	if (best.y < 40)
		best.y = (best.y + 40) / 2;
	else if (best.y > 590)
		best.y = (best.y + 590) / 2;
}

void processminis(int div, int shiftx, int shifty)
{
	nproc++;
	for(int r = MINIHW / 2 + shifty; r - MINIHW / 2 < h; r += MINIHW) {
		if (r + MINIHW / 2 >= h)
			r = h - MINIHW / 2;
		if (target == TARGET_ANCIENT) {
			if (div % 2 == 0 && r - MINIHW / 2 <= 315 && 315 <= r + MINIHW / 2 || 
				div % 3 == 0 && (r - MINIHW / 2 <= 210 && 210 <= r + MINIHW / 2 || r - MINIHW / 2 <= 420 && 420 <= r + MINIHW / 2) ||
				div % 6 == 0 && (r - MINIHW / 2 <= 105 && 105 <= r + MINIHW / 2 || r - MINIHW / 2 <= 525 && 525 <= r + MINIHW / 2))
				continue;
		} else if (div == 2 && r - MINIHW / 2 <= 315 && 315 <= r + MINIHW / 2 || 
				div == 3 && (r - MINIHW / 2 <= 210 && 210 <= r + MINIHW / 2 || r - MINIHW / 2 <= 420 && 420 <= r + MINIHW / 2))
				continue;
			for(int c = MINIHW / 2 + shiftx; c - MINIHW / 2 < w; c += MINIHW) {
				if (c + MINIHW / 2 >= w)
					c = w - MINIHW / 2;
				if (target == TARGET_ANCIENT) {
					if (div % 2 == 0 && c - MINIHW / 2 <= 618 && 618 <= c + MINIHW / 2 || 
						div % 3 == 0 && (c - MINIHW / 2 <= 412 && 412 <= c + MINIHW / 2 || c - MINIHW / 2 <= 824 && 824 <= c + MINIHW / 2) ||
						div % 6 == 0 && (c - MINIHW / 2 <= 206 && 206 <= c + MINIHW / 2 || c - MINIHW / 2 <= 1030 && 1030 <= c + MINIHW / 2))
						continue;
				} else if (div == 2 && c - MINIHW / 2 <= 618 && 618 <= c + MINIHW / 2 || 
					div == 3 && (c - MINIHW / 2 <= 412 && 412 <= c + MINIHW / 2 || c - MINIHW / 2 <= 824 && 824 <= c + MINIHW / 2))
					continue;
			pickmini(c, r);
			string s = "mini_" + name + "_bw_" + itos(div) + '_' + itos(shiftx) + '_' + itos(shifty) + '_' + itos(c) + '_' + itos(r) + ".tga";
//			savemini(s);
			hit_t best;
			if (checkall(best)) {
				writebest(best);
				char rt = resulttype(best);
				if (rt == 'r' || rt == 'R') {
					nroadriver++;
					int x = best.x + best.w / 2 + c - MINIHW / 2;
					int y = best.y + best.h / 2 + r - MINIHW / 2;
					FOR0(i, nhit) {
						if (SQR(hits[i].x - x) + SQR(hits[i].y - y) < SQR(120))
							hits[i].nroadriver++;
					}
				}
				else if (rt == 'a')
					nancient++;
				else if (rt == 'M')
					nmodern++; 
				if (target == TARGET_ANCIENT && rt == 'a' || target == TARGET_MODERN && rt == 'M') {
					best.prob = hitstats[best.type].scorerate * hitstats[best.type].goodrate;
					if (target == TARGET_ANCIENT) {
						best.prob *= (1 + min(1, 3 * (0.55 - best.err)));
						best.prob *= (1 + min(0.25, (double) max(best.h, best.w) / MINIHW / 4));
					}
					best.x += best.w / 2 + c - MINIHW / 2;
					best.y += best.h / 2 + r - MINIHW / 2;
					adjusttoedge(best);
					best.nroadriver = 0;
					hits[nhit++] = best;
					savemini(s);
				}
			}
		}
	}
}

void mergehits()
{
	if (!nhit)
		return;
//	FOR0(i, nhit)
//		cerr << hits[i].x << ' ' << hits[i].y << endl;
	// can be a much wiser partitioning, but atm it is not the big thing
	int nrealhit = 1;
	hits[0].mergecnt = 1;
	FOR(i, 1, nhit - 1) {
		bool del = false;
		FOR0(j, nrealhit) {
			int d2 = SQR(hits[i].x - hits[j].x) + SQR(hits[i].y - hits[j].y);
			if (d2 <= SQR(80)) {
				hits[j].prob = max(hits[i].prob, hits[j].prob);
				hits[j].nroadriver = max(hits[i].nroadriver, hits[j].nroadriver);
				int testx = (hits[i].x + hits[j].x * hits[j].mergecnt) / (hits[j].mergecnt + 1); // range should be based on size (and maybe prob as well)
				int testy = (hits[i].y + hits[j].y * hits[j].mergecnt) / (hits[j].mergecnt + 1);
				bool ok = true;
				FOR0(k, nrealhit) {
					if (j != k && SQR(testx - hits[k].x) + SQR(testy - hits[k].y) <=  2 * SQR(80)) {
						ok = false;
						break;
					}
				}
				if (ok) {
					hits[j].x = testx;
					hits[j].y = testy;
					hits[j].mergecnt++;
				}
				del = true;
				break;
			}
		}
		if (!del) {
			hits[i].mergecnt = 1;
			hits[nrealhit++] = hits[i];
		}
	}
	nhit = nrealhit;
}

VS assesshits()
{
	mergehits();
	FOR0(i, nhit) {
		if (hits[i].x < 300 || hits[i].y < 50 || hits[i].x > 1100 || hits[i].y > 580)
			hits[i].prob *= 0.9;
		if (target == TARGET_ANCIENT && hits[i].nroadriver >= nproc)
			hits[i].prob *= 0.66;
	}
	int nfreeslot = max(1, 5 - min(2, nroadriver / 3 / nproc) - (target == TARGET_ANCIENT ? nmodern / nproc : min(1, nancient / 3 / nproc)));
	double slotrate = min(1, (double) nfreeslot / nhit);
	if (target == TARGET_ANCIENT) {
		slotrate = min(1, (double) (3 * nfreeslot) / nhit);
	}
	cerr << "slotcalc: " << nroadriver << ' ' << nancient << ' ' << nmodern << " nproc " << nproc << " hits " << nhit << " slots " << nfreeslot << ' ' << slotrate << endl;
	VS vs;
	FOR0(i, nhit) {
		string s = itos(hits[i].x) + ',' + itos(hits[i].y) + ',' + dbltos(slotrate * hits[i].prob);
		cerr << s << ':' << hitstrs[hits[i].type] << endl;
		vs.PB(s);
	}
	return vs;
}

void initprocess()
{
	nroadriver = 0;
	nancient = 0;
	nmodern = 0;
	nhit = 0;
}

VS processall()
{
	initprocess();
	preprocess();

/*	if (target == TARGET_ANCIENT) {
		process(2, true);
		processminis(2, 0, 0);
		processminis(2, MINIHW / 2, MINIHW / 2);
	} */

	process(6);
	processminis(6, 0, 0);
	processminis(6, MINIHW / 2, MINIHW / 2);
	processminis(6, MINIHW / 2, 0);
	processminis(6, 0, MINIHW / 2);
	process(2);
	processminis(2, 0, 0);
	processminis(2, MINIHW / 2, MINIHW / 2);
	if (target != TARGET_ANCIENT || test && !timewarn1) {
		processminis(2, MINIHW / 2, 0);
		processminis(2, 0, MINIHW / 2);
	}
	if (target == TARGET_ANCIENT && test && !timewarn1) {
			processminis(2, MINIHW / 4, MINIHW / 4);
		processminis(2, 3 * MINIHW / 4, 3 * MINIHW / 4);
		processminis(2, MINIHW / 4, 3 * MINIHW / 4);
		processminis(2, 3 * MINIHW / 4, MINIHW / 4);
	} 
	process(3);
	if (target != TARGET_ANCIENT || test && !timewarn1) {
		processminis(3, 0, 0);
		processminis(3, MINIHW / 2, MINIHW / 2);
	}
	processminis(3, MINIHW / 2, 0);
	processminis(3, 0, MINIHW / 2);

	return assesshits();
}

class StructureRecognition {
public:
	string init(VS& meta, int N)
	{
		cerr << "init" << endl;
		::init();
		ntrain = SZ(meta);
		maxannotreq = N;
		processmeta(meta);
#ifdef LOCAL
		ifstream ifs("mark.txt");
		if (!ifs.eof()) {
			ifs >> train;
			cerr << "mark read: " << train << endl;
			if (train) {
				train++;
				if (mode == MODE_HITSTAT) {
					FOR0(i, MAX_HITTYPE)
						ifs >> hitstats[i].good >> hitstats[i].bad >> hitstats[i].goodannot >> hitstats[i].totalannot;
				} else {
					FOR0(i, MAX_CATEGORY)
						ifs >> catstats[i].good >> catstats[i].badmodern >> catstats[i].badancient >> catstats[i].roadriver >> catstats[i].none;
				}
			}
		}
#endif

		if (!ntrain || (mode == MODE_REAL && target != TARGET_MIX))
			return "END";
		return ids[train];
	}

	string receiveImage(string id, VI& image)
	{
		LI starttime = gettimems();
		cerr << "receive image " << id << endl;
		name = "train_" + itos(train);
		nannotblock = nannotblocks[train];
		load(image);

//		process();

#ifdef LOCAL
		string fn = "example_tile" + itos(train) + ".csv";
		if (newannotf.is_open())
			newannotf.close();
		newannotf.open(fn.c_str());
		newannotf << nannotblock << endl;
#endif
		annotblock = 0;
		ndistannot = 0;
		accumulatedtime += gettimems() - starttime;
		if (annotreq++ == maxannotreq) {
			finish();
			return "END";
		} else
			return ids[train] + ",0";
	}
   
	string receiveAnnotations(string id, int annotblock2, VS& serie) 
	{
		VS vs;
		try {
			// annot2 == annot
			LI starttime = gettimems();
			cerr << "receive annotation " <<  id << ':' << annotblock2 << '/' << annotblock << ' ' << serie[0] << endl;
			
			processserie(serie);

			while (annotblock < nannotblock - 1) {
				annotblock++;
	 			if (annotreq++ == maxannotreq) {
					finish();
					return "END";
				} else
					return ids[train] + ',' + itos(annotblock);
			}

			preprocess();
			processdistannots(false);

			if (mode == MODE_REAL && target == TARGET_MIX || mode == MODE_CHECKMINI) {
				process(3);
				processdistannots(true);
			} else if (mode == MODE_REAL && target == TARGET_ANCIENT || mode == MODE_HITSTAT) {
				processall();
				statdistannots();
			}

			accumulatedtime += gettimems() - starttime;

			if (accumulatedtime > TRAIN_TIMELIMIT) {
				cerr << "train timelimit exceeded" << endl;
				return "END";
			} 

	//!!!nondynamic ancient?
	//		if ((mode != MODE_REAL || target == TARGET_MIX)) {
			if ((mode != MODE_REAL || target != TARGET_MODERN)) {
				if (++train < ntrain)
					return ids[train];
			}
			finish();
		} catch(...) {
			cerr << "exception" << endl;
		}
		return "END";
	}
	
	VS labelImage(VI& image)
	{
		VS vs;
		try {

			LI starttime = gettimems();
			cerr << "label" << test << endl;

			if (mode == MODE_REAL && target == TARGET_MIX) {
				target = TARGET_ANCIENT;
				cerr << "just force ancient with defaults!" << endl;
			}

			name = string("test_") + itos(test);

			if (accumulatedtime > TIMELIMIT1) 
				timewarn1 = true;
			if (accumulatedtime > TIMELIMIT) {
				if (!timewarn) {
					cerr << "timelimit exceeded" << endl;
					timewarn = true;
				}
			} else {
				load(image);
				vs = processall();
			}
			test++;

			accumulatedtime += gettimems() - starttime;

			//!!! modern only
	//		if (maxannotreq == 2000)
	//			return VS();
		} catch(...) {
			cerr << "exception" << endl;
		}
		return vs;
	}
};

#ifdef LOCAL

void loadannotspecs()
{
	if (mode == MODE_REAL || mode == MODE_HITSTAT)
		return; // cant cheat
	ifstream ifs("annot_spec.txt");
	while(!ifs.eof()) {
		int dummy;
		int t;
		annotspec_t spec;
		ifs >> spec.category >> t >> dummy >> dummy >> spec.x >> spec.y;
		if (spec.category <= '9')
			continue;
		annotspecs[t][nannotspecs[t]++] = spec;
	}
}

void testblock();

int main(int argc, char *argv[])
{
	srand(1);
	debugf = stderr;

	loadannotspecs();
	
	if (argc >= 2)
		mode = atoi(argv[1]);

	if (argc >= 3)
		target = atoi(argv[2]);

	if (mode == MODE_SINGLE) {
		testblock();
		return 1;
	}

	StructureRecognition sr;

	int d, size;
	string s, id;
	VS vs;
	VI vi;

    cin >> d;
    cin >> size;
    getline(cin, s);
	FOR0(i, size) {
        getline(cin, s);
		vs.PB(s);
	}

    string cmd = sr.init(vs, d);
    cout << cmd << endl;
    cout.flush();
    
    while (true) {
        int visualcmd;
		cin >> visualcmd;
		cerr << visualcmd << endl;
		getline(cin, s);
		switch (visualcmd) {
		case 0:
            return 1;
		case 1:
			getline(cin, id);
			cin >> size;
			vi.clear();
			FOR0(i, size) {
				cin >> d;
				vi.PB(d);
			}
            cmd = sr.receiveImage(id, vi);
            cout << cmd << endl;
			cout.flush();
			break;
		case 2:
			getline(cin, id);
            int annotator;
			cin >> annotator;
			cin >> size;
			getline(cin, s);
			vs.clear();
			FOR0(i, size) {
				getline(cin, s);
				vs.PB(s);
			}
            cmd = sr.receiveAnnotations(id, annotator, vs);
            cout << cmd << endl;
			cout.flush();
			break;
		case 3:
			cin >> size;
			vi.clear();
			FOR0(i, size) {
				cin >> d;
				vi.PB(d);
			}
            VS labels = sr.labelImage(vi);
            cout << SZ(labels) << endl;
			FOR0(i, SZ(labels))
				cout << labels[i] << endl;
			cout.flush();
        }
	}
}

void testblock()
{
	init();
	target = TARGET_ANCIENT;
//	target = TARGET_MODERN;
	loaddat("train_13.dat");
#if 1
	processall();
#else
	preprocess();
	process(3);
	pickmini(994, 531);
	savemini("mini_bw.tga");
	calcbest('r'); 
#endif
}

#endif

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