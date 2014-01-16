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

// TopCoder Marathon - StructureRecognition
// author: elder1g



#include <cctype>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <string>
#include <string.h>
#include <sstream>
#include <iostream>
#include <fstream>
#include <map>
#include <set>
#include <deque>
#include <sys/time.h>

#ifdef TEST_IMG
#include "img_lib.cpp"
#endif

using namespace std;

#ifdef LOCAL
	ofstream cmd_rec("cmd_rec.txt");
#else
	ostream& cmd_rec=cerr;
#endif

double until_now=0;

#define MAX_NODE_CNT 50000
#define MAX_LEVEL 15
#define TRAIN_IMG_CNT 839
#define IMG_WIDTH 1236
#define IMG_HEIGHT 630
#define IMG_BLOCK_WIDTH 1237
#define IMG_BLOCK_HEIGHT 631
#define BLUR_CNT 3
#define GREY_LEVEL_CNT 8
#define GRADIENT_KIND_CNT 6
#define ANGLE_LEVEL_CNT 8
#define RARE_COLOR_LEVEL 5
#define RARE_COLOR_KIND 16

vector<int> hold_id[MAX_NODE_CNT];
const double eps=1e-12;

int img_r[BLUR_CNT][IMG_HEIGHT][IMG_WIDTH];
int img_g[BLUR_CNT][IMG_HEIGHT][IMG_WIDTH];
int img_b[BLUR_CNT][IMG_HEIGHT][IMG_WIDTH];
int img_temp_r[IMG_HEIGHT][IMG_WIDTH];
int img_temp_g[IMG_HEIGHT][IMG_WIDTH];
int img_temp_b[IMG_HEIGHT][IMG_WIDTH];
int img_grey[BLUR_CNT][IMG_HEIGHT][IMG_WIDTH];
int tmp_grey2[IMG_HEIGHT][IMG_WIDTH];
int img_grey_block[BLUR_CNT][IMG_BLOCK_HEIGHT][IMG_BLOCK_WIDTH];
long long img_grey2_block[BLUR_CNT][IMG_BLOCK_HEIGHT][IMG_BLOCK_WIDTH];
int img_tot_color[BLUR_CNT][4];
int img_grey_level[GREY_LEVEL_CNT][IMG_HEIGHT][IMG_WIDTH];
int img_grey_level_block[GREY_LEVEL_CNT][IMG_BLOCK_HEIGHT][IMG_BLOCK_WIDTH];
int img_angle_level[ANGLE_LEVEL_CNT][IMG_HEIGHT][IMG_WIDTH];
int img_angle_level_block[ANGLE_LEVEL_CNT][IMG_BLOCK_HEIGHT][IMG_BLOCK_WIDTH];
/*
int tmp_color_group[IMG_HEIGHT][IMG_WIDTH];
int img_rare_level[RARE_COLOR_LEVEL][IMG_HEIGHT][IMG_WIDTH];
int img_rare_level_block[RARE_COLOR_LEVEL][IMG_BLOCK_HEIGHT][IMG_BLOCK_WIDTH];
double nowr[RARE_COLOR_KIND];
double nowg[RARE_COLOR_KIND];
double nowb[RARE_COLOR_KIND];
*/

/*
int img_red_level[BLUR_CNT][GREY_LEVEL_CNT][IMG_HEIGHT][IMG_WIDTH];
int img_green_level[BLUR_CNT][GREY_LEVEL_CNT][IMG_HEIGHT][IMG_WIDTH];
int img_blue_level[BLUR_CNT][GREY_LEVEL_CNT][IMG_HEIGHT][IMG_WIDTH];
long long img_red_level_block[BLUR_CNT][GREY_LEVEL_CNT][IMG_BLOCK_HEIGHT][IMG_BLOCK_WIDTH];
long long img_green_level_block[BLUR_CNT][GREY_LEVEL_CNT][IMG_BLOCK_HEIGHT][IMG_BLOCK_WIDTH];
long long img_blue_level_block[BLUR_CNT][GREY_LEVEL_CNT][IMG_BLOCK_HEIGHT][IMG_BLOCK_WIDTH];
*/

int img_gradient_x[IMG_HEIGHT][IMG_WIDTH];
int img_gradient_y[IMG_HEIGHT][IMG_WIDTH];
int img_gradient_all[IMG_HEIGHT][IMG_WIDTH];
int img_gradient_level[GRADIENT_KIND_CNT][IMG_HEIGHT][IMG_WIDTH];
int img_gradient_level_block[GRADIENT_KIND_CNT][IMG_BLOCK_HEIGHT][IMG_BLOCK_WIDTH];
double PI=acos(-1.0);

int maybe_cnt[500];
int select_cnt[500];

int offr8[8]={-1,-1,-1,0,0,1,1,1};
int offc8[8]={-1,0,1,-1,1,-1,0,1};
struct GraGroup
{
	int cnt;
	int maxx,minx,maxy,miny;
};
int now_group_index=0;
vector<GraGroup> gra_groups;
int img_gradient_indexes[IMG_HEIGHT][IMG_WIDTH];
/*
int r_asn[100000];
int g_asn[100000];
int b_asn[100000];
*/

bool asn_wts=false;
vector<int> forest_wts;

string convert_int(int a)
{
	stringstream sa;
	sa<<a;
	return sa.str();
}

string convert_double(double a)
{
	stringstream sa;
	sa<<a;
	return sa.str();
}

int read_int(string a)
{
	stringstream sa(a);
	int t;
	sa>>t;
	return t;
}

struct Annotator
{
	int ord;
	int before;
	int all;
	int seconds;
	Annotator()
	{
	}
	Annotator(int od,int bef,int al,int sec)
	{
		ord=od;
		before=bef;
		all=al;
		seconds=sec;
	}
};

bool operator < (const Annotator a,const Annotator b)
{
	return a.before>b.before;
}

struct AntRec
{
	int x;
	int y;
	int ord;
	AntRec()
	{
	}
	AntRec(int _x,int _y,int _ord)
	{
		x=_x;
		y=_y;
		ord=_ord;
	}
};

int N;
vector<Annotator> annotators[TRAIN_IMG_CNT];
string train_img_name[TRAIN_IMG_CNT];
vector<pair<int, int> > chk_train_img_order;
int to_chk_img_idx;
int each_image_chk_cnt;
int each_image_collect_cnt;
int left_aloc;
int to_chk_img_ord;
int to_check_ant_idx;
vector<AntRec> antrec;
vector<double> fact_pool;
vector<vector<int> > prop_pool;
int last_mark_ord[1000000];
int atr_idx;
bool build_label=false;
int label_index=0;

struct TreeNode
{
	int left; //right=left+1
	int split_index; //leaf -> -1
	double split_value; //leaf -> prediction
	int level;
};

unsigned int m1=11;
unsigned int m2=102;

void resetrand()
{
	m1=11;
	m2=102;
}
unsigned int quickrand()
{
    m1=36969*(m1&65535)+(m1>>16);
    m2=18000*(m2&65535)+(m2>>16);
    return (m1<<16)+m2;
}

struct BinaryTree
{
	vector<TreeNode> nds;
	double make_decision(const vector<int>& properties)
	{
		int cur_node=0;
		while(nds[cur_node].split_index>=0)
		{
			double v0=nds[cur_node].split_value;
			int v1=properties[nds[cur_node].split_index];
			if (v1<v0) cur_node=nds[cur_node].left;
			else cur_node=nds[cur_node].left+1;
		}
		return nds[cur_node].split_value;
	}
	double get_prd(const vector<double>& fact,const vector<int>& ids)
	{
		double sum=0;
		for (int i=0;i<ids.size();++i)
			sum+=fact[ids[i]];
		return sum/ids.size();
	}
	void construct_tree(const vector<double>& fact,const vector<vector<int> >& properties)
	{
		double fact_mark=0;
		unsigned int prop_mark=0;
		for(int i=0;i<fact.size();++i)
			fact_mark+=fact[i];
		for(int i=0;i<properties.size();++i)
			for(int j=0;j<properties[i].size();++j)
				prop_mark+=properties[i][j];
		cerr<<"Mark all:"<<fact_mark<<" "<<prop_mark<<endl;			
		int todeal=0;
		nds.reserve(MAX_NODE_CNT);
		nds.push_back(TreeNode());
		int row_size=fact.size();
		int col_size=properties[0].size();
		int sel_size=6;
		TreeNode& root=nds[0];
		fact_mark=0;
		prop_mark=0;
		vector<int> cpool;
		for (int i=0;i<col_size;++i)
			for (int j=0;j<forest_wts[i];++j)
				cpool.push_back(i);
		for (int i=0;i<row_size;++i)//try smaller sample size(than pool size?)
		{
			int ns=quickrand()%row_size;
			//int ns=i;
			hold_id[0].push_back(ns);
			fact_mark+=fact[ns];
			for (int j=0;j<properties[ns].size();++j)
				prop_mark+=properties[ns][j];
		}
		cerr<<"Mark choose:"<<fact_mark<<" "<<prop_mark<<endl;	
		root.level=1;
		/*
		for (int i=0;i<hold_id[0].size();++i)
		{
			int ci=hold_id[0][i];
			if(fact[ci]==0) continue;
			cerr<<fact[ci]<<" ";
			for(int j=0;j<24;++j)
				cerr<<properties[ci][j]<<" ";
			cerr<<endl;
		}
		*/
		int end_split=100;
		while(todeal<nds.size())
		{
			vector<int> chs_from;
			TreeNode& node=nds[todeal];
			bool leaf=false;
			if (hold_id[todeal].size()<end_split || node.level==MAX_LEVEL)
				leaf=true;
			else
			{
				chs_from.reserve(sel_size);
				for (int i=0;i<sel_size;++i)
				{
					int ns=quickrand()%cpool.size();
					chs_from.push_back(cpool[ns]);
					maybe_cnt[cpool[ns]]++;
				}
				double min_var=1e+100;
				int min_id=-1;
				double min_sp_value=-1;
				int best_left_cnt=-1;
				int best_right_cnt=-1;
				for (int i=0;i<chs_from.size();++i)
				{
					vector<pair<int,double> > sortlist;
					sortlist.reserve(hold_id[todeal].size());
					int pro_id=chs_from[i];
					double left_sumx=0;
					double right_sumx=0;
					double left_sumx2=0;
					double right_sumx2=0;
					int left_cnt=0;
					int right_cnt=hold_id[todeal].size();
					for (int j=0;j<hold_id[todeal].size();++j)
					{
						int idx=hold_id[todeal][j];
						double vle=fact[idx];
						sortlist.push_back(make_pair(properties[idx][pro_id],vle));
						right_sumx+=vle;
						right_sumx2+=vle*vle;
					}
					sort(sortlist.begin(),sortlist.end());
										
					for (int j=0;j+1<sortlist.size();++j)
					{
						int v1=sortlist[j].second;
						//int v2=sortlist[j+1].second;
						left_sumx+=v1;
						right_sumx-=v1;
						left_sumx2+=v1*v1;
						right_sumx2-=v1*v1;
						++left_cnt;
						--right_cnt;
						//if (fabs(v1-v2)<eps) continue;
						if (sortlist[j].first==sortlist[j+1].first) continue;
						if (left_cnt<end_split || right_cnt<end_split) continue;
						double left_ex=left_sumx/left_cnt;
						double left_ex2=left_sumx2/left_cnt;
						double left_var=left_ex2-left_ex*left_ex;
						
						if(left_var<0)
							left_var=0;
						else
							left_var=sqrt(left_var*left_cnt/(left_cnt-1));
						
						double right_ex=right_sumx/right_cnt;
						double right_ex2=right_sumx2/right_cnt;
						double right_var=right_ex2-right_ex*right_ex;
						
						if(right_var<0)
							right_var=0;
						else
							right_var=sqrt(right_var*right_cnt/(right_cnt-1));						
						
						double ova_var=left_var*left_cnt+right_var*right_cnt;
						
						if (ova_var*(1+eps)<min_var)
						{
							min_sp_value=(sortlist[j].first+sortlist[j+1].first)*0.5;
							min_var=ova_var;
							min_id=pro_id;
							best_left_cnt=left_cnt;
							best_right_cnt=right_cnt;
						}
					}
				}
				if (min_var<1e+80)
				{
					select_cnt[min_id]++;
					node.split_index=min_id;
					node.split_value=min_sp_value;
					node.left=nds.size();
					TreeNode tleft,tright;
					nds.push_back(tleft);
					nds.push_back(tright);
					nds[node.left].level=nds[node.left+1].level=node.level+1;
					hold_id[node.left].reserve(best_left_cnt);
					hold_id[node.left+1].reserve(best_right_cnt);
					for (int i=0;i<hold_id[todeal].size();++i)
						if (properties[hold_id[todeal][i]][node.split_index]<node.split_value)
							hold_id[node.left].push_back(hold_id[todeal][i]);
						else
							hold_id[node.left+1].push_back(hold_id[todeal][i]);
				}
				else
					leaf=true;
			}
			if (leaf)
			{
				node.split_index=-1;
				node.split_value=get_prd(fact,hold_id[todeal]);
			}
			++todeal;
		}
		//clear data
		for (int i=0;i<MAX_NODE_CNT;++i)
		{
			vector<int> empt;
			hold_id[i].swap(empt);
		}
	}
};

#define TREES_CNT 40
BinaryTree trees[TREES_CNT];

void build_trees()
{
	for (int i=0;i<TREES_CNT;++i)
	{
		m1=11+i;
		m2=102+i;
		trees[i].construct_tree(fact_pool,prop_pool);
		cerr<<fact_pool.size()<<" "<<prop_pool.size()<<" ";
		cerr<<"Build "<<i<<":"<<trees[i].nds.size()<<endl;
		/*
		for(int j=0;j<trees[i].nds.size();++j)
			cerr<<trees[i].nds[j].split_index<<" ";
		cerr<<endl;
		*/
	}
}

struct StructureRecognition
{
	string init(vector<string>& imageMetaData,int _N)
	{
		N=_N;
		resetrand();
		/*
		for (int i=0;i<100000;++i)
		{
			r_asn[i]=quickrand()%256;
			g_asn[i]=quickrand()%256;
			b_asn[i]=quickrand()%256;
		}
		*/
		resetrand();
		cerr<<"Allow request:"<<N<<endl;
		//vector<Annotator> annotators[TRAIN_IMG_CNT];
		//string train_img_name[TRAIN_IMG_CNT];
		for (int i=0;i<imageMetaData.size();++i)
		{
			string str=imageMetaData[i];
			for (int j=0;j<str.size();++j)
				if (str[j]==',')
					str[j]=' ';
			stringstream sa(str);
			sa>>train_img_name[i];
			annotators[i].reserve(16);
			Annotator ant;
			int ant_ord=0;
			while(sa>>ant.before>>ant.all>>ant.seconds)
			{
				ant.ord=ant_ord;
				annotators[i].push_back(ant);
				ant_ord++;
			}
			
			for (int j=0;j<annotators[i].size();++j)
			{
				int sf=j+quickrand()%(annotators[i].size()-j);
				swap(annotators[i][j],annotators[i][sf]);
			}
			
			//sort(annotators[i].begin(),annotators[i].end());
		}
		//may add: sort each images' annotators based on past history
		resetrand();
		
		chk_train_img_order.reserve(imageMetaData.size());
		for (int i=0;i<imageMetaData.size();++i)
			chk_train_img_order.push_back(make_pair(-annotators[i].size(),i));
		sort(chk_train_img_order.begin(),chk_train_img_order.end());
		
		for (int i=0;i<imageMetaData.size();++i)
		{
			int sf=i+quickrand()%(imageMetaData.size()-i);
			swap(chk_train_img_order[i],chk_train_img_order[sf]);
		}
		resetrand();
		
		to_chk_img_ord=0;
		to_chk_img_idx=chk_train_img_order[to_chk_img_ord].second;
		/*
		each_image_chk_cnt=(int)((N/TRAIN_IMG_CNT)+1);
		each_image_chk_cnt=min(each_image_chk_cnt,10);
		each_image_chk_cnt=max(each_image_chk_cnt,5);
		*/
		each_image_chk_cnt=10;//5
		//if(N>=2000)
		//	each_image_chk_cnt=10;
		each_image_collect_cnt=1400000/(N/each_image_chk_cnt);
		left_aloc=N;
		fact_pool.reserve(2000000);
		prop_pool.reserve(2000000);
		return train_img_name[to_chk_img_idx];
	}
	void generate_block(int raw[IMG_HEIGHT][IMG_WIDTH],long long block[IMG_BLOCK_HEIGHT][IMG_BLOCK_WIDTH],int szof)
	{
		memset(block,0,szof);
		for (int r=1;r<IMG_BLOCK_HEIGHT;++r)
		{
			int tl=0;
			for(int c=1;c<IMG_BLOCK_WIDTH;++c)
			{
				tl+=raw[r-1][c-1];
				block[r][c]=block[r-1][c]+tl;
			}
		}
	}
	void generate_block_int_ver(int raw[IMG_HEIGHT][IMG_WIDTH],int block[IMG_BLOCK_HEIGHT][IMG_BLOCK_WIDTH],int szof)
	{
		memset(block,0,szof);
		for (int r=1;r<IMG_BLOCK_HEIGHT;++r)
		{
			int tl=0;
			for(int c=1;c<IMG_BLOCK_WIDTH;++c)
			{
				tl+=raw[r-1][c-1];
				block[r][c]=block[r-1][c]+tl;
			}
		}
	}
	void generate_image_info(vector<int>& image)
	{
		int idx=2;
		for (int r=0;r<IMG_HEIGHT;++r)
		{
			for (int c=0;c<IMG_WIDTH;++c)
			{
				img_r[0][r][c]=(image[idx]>>16);
				img_g[0][r][c]=((image[idx]>>8)&255);
				img_b[0][r][c]=(image[idx]&255);
				idx++;
			}
		}
		memcpy(img_temp_r,img_r[0],sizeof(img_temp_r));
		memcpy(img_temp_g,img_g[0],sizeof(img_temp_g));
		memcpy(img_temp_b,img_b[0],sizeof(img_temp_b));
		for (int i=0;i<3;++i)
		{
			for (int r=0;r<IMG_HEIGHT;++r)
			{
				for (int c=0;c<IMG_WIDTH;++c)
				{
					double ct=0;
					double sum_r=0;
					double sum_g=0;
					double sum_b=0;
					for (int dr=-1;dr<=1;++dr)
					{
						for (int dc=-1;dc<=1;++dc)
						{
							int nr=r+dr;
							int nc=c+dc;
							if (nr>=0 && nr<IMG_HEIGHT && nc>=0 && nc<IMG_WIDTH)
							{
								sum_r+=img_temp_r[nr][nc];
								sum_g+=img_temp_g[nr][nc];
								sum_b+=img_temp_b[nr][nc];
								++ct;
							}
						}
					}
					int tr=(int)(sum_r/ct+0.5);
					tr=min(255,tr);
					int tg=(int)(sum_g/ct+0.5);
					tg=min(255,tg);
					int tb=(int)(sum_b/ct+0.5);
					tb=min(255,tb);
					img_r[1][r][c]=tr;
					img_g[1][r][c]=tg;
					img_b[1][r][c]=tb;
				}
			}
			memcpy(img_temp_r,img_r[1],sizeof(img_temp_r));
			memcpy(img_temp_g,img_g[1],sizeof(img_temp_g));
			memcpy(img_temp_b,img_b[1],sizeof(img_temp_b));
		}
		//memcpy(img_r[2],img_r[1],sizeof(img_r[2]));
		//memcpy(img_g[2],img_g[1],sizeof(img_g[2]));
		//memcpy(img_b[2],img_b[1],sizeof(img_b[2]));
		//60,70,115->60,70,100->50,70,100->50,60,100
		int threshold_dir1a=50;
		//int threshold_dir1b=70;
		int threshold_dir2a=60;
		int threshold_dir2b=100;
		//int grey_cnt[256];
		//int map_grey[256];
		for (int m=0;m<2;++m)
		{
			if (m==0)
				memset(img_grey_level,0,sizeof(img_grey_level));
			/*
			memset(img_red_level[m],0,sizeof(img_red_level[m]));
			memset(img_green_level[m],0,sizeof(img_green_level[m]));
			memset(img_blue_level[m],0,sizeof(img_blue_level[m]));
			*/
			if (m==1)
				memset(img_gradient_level,0,sizeof(img_gradient_level));
			for (int r=0;r<IMG_HEIGHT;++r)
			{
				for (int c=0;c<IMG_WIDTH;++c)
				{
					img_grey[m][r][c]=(int)((img_r[m][r][c]*299+img_g[m][r][c]*587+img_b[m][r][c]*114+500)/1000.0);
					img_grey[m][r][c]=min(255,img_grey[m][r][c]);
				}
			}
			/*
			if (m==2)
			{
				memset(grey_cnt,0,sizeof(grey_cnt));
				int below_cnt=0;
				for (int r=0;r<IMG_HEIGHT;++r)
					for (int c=0;c<IMG_WIDTH;++c)
						grey_cnt[img_grey[m][r][c]]++;
				double tot_point_cnt=IMG_HEIGHT*IMG_WIDTH;
				for (int i=0;i<256;++i)
				{
					int place=(int)((below_cnt+grey_cnt[i]/2.0)/tot_point_cnt*255+0.5);
					place=min(place,255);
					map_grey[i]=place;
					below_cnt+=grey_cnt[i];
				}
				for (int r=0;r<IMG_HEIGHT;++r)
					for (int c=0;c<IMG_WIDTH;++c)
						img_grey[m][r][c]=map_grey[img_grey[m][r][c]];
			}
			*/
			memset(img_tot_color[m],0,sizeof(img_tot_color[m]));
			for (int r=0;r<IMG_HEIGHT;++r)
				for (int c=0;c<IMG_WIDTH;++c)
				{
					img_tot_color[m][0]+=img_grey[m][r][c];
					img_tot_color[m][1]+=img_r[m][r][c];
					img_tot_color[m][2]+=img_g[m][r][c];
					img_tot_color[m][3]+=img_b[m][r][c];
				}
			if (m==0)
			{
				for (int r=0;r<IMG_HEIGHT;++r)
				{
					for (int c=0;c<IMG_WIDTH;++c)
					{
						int glv=img_grey[m][r][c]/(256/GREY_LEVEL_CNT);
						img_grey_level[glv][r][c]++;
						/*
						glv=img_r[m][r][c]/(256/GREY_LEVEL_CNT);
						img_red_level[m][glv][r][c]++;
						glv=img_g[m][r][c]/(256/GREY_LEVEL_CNT);
						img_green_level[m][glv][r][c]++;
						glv=img_b[m][r][c]/(256/GREY_LEVEL_CNT);
						img_blue_level[m][glv][r][c]++;
						*/
					}
				}
			
				for (int i=0;i<GREY_LEVEL_CNT;++i)
				{
					generate_block_int_ver(img_grey_level[i],img_grey_level_block[i],sizeof(img_grey_level_block[i]));
					/*
					generate_block_int_ver(img_red_level[m][i],img_red_level_block[m][i],sizeof(img_red_level_block[m][i]));
					generate_block_int_ver(img_green_level[m][i],img_green_level_block[m][i],sizeof(img_green_level_block[m][i]));
					generate_block_int_ver(img_blue_level[m][i],img_blue_level_block[m][i],sizeof(img_blue_level_block[m][i]));
					*/
				}
			}
			
			generate_block_int_ver(img_grey[m],img_grey_block[m],sizeof(img_grey_block[m]));
			for (int r=0;r<IMG_HEIGHT;++r)
				for (int c=0;c<IMG_WIDTH;++c)
					tmp_grey2[r][c]=img_grey[m][r][c]*img_grey[m][r][c];
			generate_block(tmp_grey2,img_grey2_block[m],sizeof(img_grey2_block[m]));
			
			//int img_gradient_x[IMG_HEIGHT][IMG_WIDTH];
			//int img_gradient_y[IMG_HEIGHT][IMG_WIDTH];
			//int img_gradient_all[IMG_HEIGHT][IMG_WIDTH];
			if (m==1)
			{
				for (int r=0;r<IMG_HEIGHT;++r)
				{
					for (int c=0;c<IMG_WIDTH;++c)
					{
						img_gradient_x[r][c]=img_gradient_y[r][c]=0;
						if (r>=1 && r+1<IMG_HEIGHT && c>=1 && c+1<IMG_WIDTH)
						{
							img_gradient_x[r][c]=img_grey[m][r][c+1]*2+img_grey[m][r+1][c+1]+img_grey[m][r-1][c+1]-img_grey[m][r][c-1]*2-img_grey[m][r+1][c-1]-img_grey[m][r-1][c-1];
							img_gradient_y[r][c]=img_grey[m][r+1][c]*2+img_grey[m][r+1][c+1]+img_grey[m][r+1][c-1]-img_grey[m][r-1][c]*2-img_grey[m][r-1][c-1]-img_grey[m][r-1][c+1];
						}
						img_gradient_all[r][c]=(int)(sqrt(img_gradient_x[r][c]*img_gradient_x[r][c]*1.0+img_gradient_y[r][c]*img_gradient_y[r][c])+0.5);
						if(img_gradient_x[r][c]>threshold_dir1a)
							img_gradient_level[0][r][c]++;
						if(img_gradient_x[r][c]<-threshold_dir1a)
							img_gradient_level[1][r][c]++;
						if(img_gradient_y[r][c]>threshold_dir1a)
							img_gradient_level[2][r][c]++;
						if(img_gradient_y[r][c]<-threshold_dir1a)
							img_gradient_level[3][r][c]++;
						if(img_gradient_all[r][c]>threshold_dir2a)
							img_gradient_level[4][r][c]++;
						if(img_gradient_all[r][c]>threshold_dir2b)
							img_gradient_level[5][r][c]++;
					}
				}
				//img_gradient_level[GRADIENT_KIND_CNT][IMG_HEIGHT][IMG_WIDTH];
				for (int i=0;i<GRADIENT_KIND_CNT;++i)
					generate_block_int_ver(img_gradient_level[i],img_gradient_level_block[i],sizeof(img_gradient_level_block[i]));
					
				//int img_angle_level[ANGLE_LEVEL_CNT][IMG_HEIGHT][IMG_WIDTH];
				//int img_angle_level_block[ANGLE_LEVEL_CNT][IMG_BLOCK_HEIGHT][IMG_BLOCK_WIDTH];
				memset(img_angle_level,0,sizeof(img_angle_level));
				for (int r=0;r<IMG_HEIGHT;++r)
				{
					for (int c=0;c<IMG_WIDTH;++c)
					{
						if (img_gradient_all[r][c]>threshold_dir2a)
						{
							double ang=atan2(img_gradient_y[r][c],img_gradient_x[r][c]);
							while(ang<0) ang+=2*PI;
							while(ang>2*PI) ang-=2*PI;
							int lv=ang/(2*PI/ANGLE_LEVEL_CNT);
							lv=min(lv,ANGLE_LEVEL_CNT-1);
							lv=max(lv,0);
							img_angle_level[lv][r][c]++;
						}
					}
				}
				for (int i=0;i<ANGLE_LEVEL_CNT;++i)
					generate_block_int_ver(img_angle_level[i],img_angle_level_block[i],sizeof(img_angle_level_block[i]));
				//#define RARE_COLOR_KIND 10
				//int img_rare_level[RARE_COLOR_LEVEL][IMG_HEIGHT][IMG_WIDTH];
				//int img_rare_level_block[RARE_COLOR_LEVEL][IMG_BLOCK_HEIGHT][IMG_BLOCK_WIDTH];
				//int tmp_color_group[IMG_HEIGHT][IMG_WIDTH];
				
				/*
				double totr[RARE_COLOR_KIND];
				double totg[RARE_COLOR_KIND];
				double totb[RARE_COLOR_KIND];
				int cnt_each_group[RARE_COLOR_KIND];
				double accu_group[RARE_COLOR_KIND];
				for (int i=0;i<RARE_COLOR_KIND;++i)
				{
					bool rep=false;
					do
					{
						int nr=quickrand()%IMG_HEIGHT;
						int nc=quickrand()%IMG_WIDTH;
						nowr[i]=img_r[1][nr][nc];
						nowg[i]=img_g[1][nr][nc];
						nowb[i]=img_b[1][nr][nc];
						rep=false;
						for (int j=0;j<i;++j)
							if (nowr[i]==nowr[j] && nowg[i]==nowg[j] && nowb[i]==nowb[j])
							{
								rep=true;
								break;
							}
					}while(rep);
				}
				for (int itr=0;itr<7;++itr)
				{
					for (int r=0;r<IMG_HEIGHT;++r)
					{
						for (int c=0;c<IMG_WIDTH;++c)
						{
							double mindis=1e+300;
							int gid=-1;
							for (int g=0;g<RARE_COLOR_KIND;++g)
							{
								double disr=nowr[g]-img_r[1][r][c];
								double disg=nowg[g]-img_g[1][r][c];
								double disb=nowb[g]-img_b[1][r][c];
								double tdis=disr*disr+disg*disg+disb*disb;
								if (tdis*(1+eps)<mindis)
								{
									mindis=tdis;
									gid=g;
								}
							}
							tmp_color_group[r][c]=gid;
						}
					}
					for (int i=0;i<RARE_COLOR_KIND;++i)
					{
						memset(totr,0,sizeof(totr));
						memset(totg,0,sizeof(totg));
						memset(totb,0,sizeof(totb));
						memset(cnt_each_group,0,sizeof(cnt_each_group));
					}
					for (int r=0;r<IMG_HEIGHT;++r)
					{
						for (int c=0;c<IMG_WIDTH;++c)
						{
							int gid=tmp_color_group[r][c];
							totr[gid]+=img_r[1][r][c];
							totg[gid]+=img_g[1][r][c];
							totb[gid]+=img_b[1][r][c];
							cnt_each_group[gid]++;
						}
					}
					for (int i=0;i<RARE_COLOR_KIND;++i)
					{
						if (cnt_each_group[i]==0)
						{
							nowr[i]=quickrand()%256;
							nowg[i]=quickrand()%256;
							nowb[i]=quickrand()%256;
						}
						else
						{
							nowr[i]=totr[i]/cnt_each_group[i];
							nowg[i]=totg[i]/cnt_each_group[i];
							nowb[i]=totb[i]/cnt_each_group[i];
						}
					}
				}
				vector<pair<double,int> > sort_color;
				for (int i=0;i<RARE_COLOR_KIND;++i)
					sort_color.push_back(make_pair(cnt_each_group[i]*1.0/IMG_HEIGHT/IMG_WIDTH,i));
				sort(sort_color.begin(),sort_color.end());
				double now_accu=0;
				for (int i=0;i<RARE_COLOR_KIND;++i)
				{
					now_accu+=sort_color[i].first;
					cerr<<now_accu<<" ";
					accu_group[sort_color[i].second]=now_accu;
				}
				cerr<<endl;
				memset(img_rare_level,0,sizeof(img_rare_level));
				for (int r=0;r<IMG_HEIGHT;++r)
				{
					for (int c=0;c<IMG_WIDTH;++c)
					{
						double accu=accu_group[tmp_color_group[r][c]];
						if (accu<0.05) img_rare_level[0][r][c]++;
						if (accu<0.10) img_rare_level[1][r][c]++;
						if (accu<0.20) img_rare_level[2][r][c]++;
						if (accu<0.30) img_rare_level[3][r][c]++;
						if (accu<0.50) img_rare_level[4][r][c]++;
					}
				}
				for (int i=0;i<RARE_COLOR_LEVEL;++i)
					generate_block_int_ver(img_rare_level[i],img_rare_level_block[i],sizeof(img_rare_level_block[i]));
				*/
			}
		}
		
		now_group_index=0;
		gra_groups.clear();
		memset(img_gradient_indexes,-1,sizeof(img_gradient_indexes));
		for (int r=0;r<IMG_HEIGHT;++r)
		{
			for (int c=0;c<IMG_WIDTH;++c)
			{ 
				if (img_gradient_indexes[r][c]<0 && img_gradient_all[r][c]>threshold_dir2a)
				{
					deque<pair<int,int> > q;
					q.push_back(make_pair(r,c));
					img_gradient_indexes[r][c]=now_group_index;
					GraGroup gp;
					gp.cnt=0;
					gp.maxx=-1;
					gp.maxy=-1;
					gp.minx=2000;
					gp.miny=2000;
					while(!q.empty())
					{
						int cr=q.front().first;
						int cc=q.front().second;
						gp.cnt++;
						gp.maxx=max(gp.maxx,cc);
						gp.minx=min(gp.minx,cc);
						gp.maxy=max(gp.maxy,cr);
						gp.miny=min(gp.miny,cr);
						q.pop_front();
						for (int k=0;k<8;++k)
						{
							int nr=cr+offr8[k];
							int nc=cc+offc8[k];
							if (nr>=0 && nc>=0 && nr<IMG_HEIGHT && nc<IMG_WIDTH && img_gradient_indexes[nr][nc]<0 && img_gradient_all[nr][nc]>threshold_dir2a)
							{
								img_gradient_indexes[nr][nc]=now_group_index;
								q.push_back(make_pair(nr,nc));
							}
						}
					}
					gra_groups.push_back(gp);
					now_group_index++;
				}
			}
		}
		/*
		cerr<<gra_groups.size()<<endl;
		int tid=-1;
		int tmax=-1;
		for (int i=0;i<gra_groups.size();++i)
		{
			if (gra_groups[i].cnt>tmax)
			{
				tmax=gra_groups[i].cnt;
				tid=i;
			}
		}
		if (tid>=0)
			cerr<<gra_groups[tid].cnt<<" "<<gra_groups[tid].minx<<" "<<gra_groups[tid].maxx<<" "<<gra_groups[tid].miny<<" "<<gra_groups[tid].maxy<<endl;
		*/
	}
	int get_block_value(long long block[IMG_BLOCK_HEIGHT][IMG_BLOCK_WIDTH],int r1,int c1,int r2,int c2)
	{
		long long a=block[r2+1][c2+1];
		long long b=block[r1][c2+1];
		long long c=block[r2+1][c1];
		long long d=block[r1][c1];
		return (int)(a-b-c+d);
	}
	int get_block_value_int_ver(int block[IMG_BLOCK_HEIGHT][IMG_BLOCK_WIDTH],int r1,int c1,int r2,int c2)
	{
		long long a=block[r2+1][c2+1];
		long long b=block[r1][c2+1];
		long long c=block[r2+1][c1];
		long long d=block[r1][c1];
		return (int)(a-b-c+d);
	}
	vector<int> get_prop(int r,int c)
	{
		vector<int> ret;
		if (prop_pool.size()!=0)
			ret.reserve(prop_pool[0].size());
		/*
		int rg0[4]={10,20,30,40};
		for (int i=1;i<4;++i)
		{
			int j=i-1;
			{
				int out_tx=get_block_value_int_ver(img_grey_block[0],r-rg0[i],c-rg0[i],r+rg0[i],c+rg0[i]);
				int inn_tx=get_block_value_int_ver(img_grey_block[0],r-rg0[j],c-rg0[j],r+rg0[j],c+rg0[j]);
				int points_out=rg0[i]*2+1;
				int points_inn=rg0[j]*2+1;
				double ex1=1.0*(out_tx-inn_tx)/(points_out-points_inn);
				double ex2=1.0*inn_tx/points_inn;
				ret.push_back((int)((ex1-ex2)*100));
				//double ex=1.0*out_tx/points_out;
				//ret.push_back((int)(ex*100));
				if (!asn_wts) forest_wts.push_back(1);
			}
		}
		*/
		//int img_angle_level[ANGLE_LEVEL_CNT][IMG_HEIGHT][IMG_WIDTH];
		//int img_angle_level_block[ANGLE_LEVEL_CNT][IMG_BLOCK_HEIGHT][IMG_BLOCK_WIDTH];
		/*
		int rga[4]={10,20,30,40};
		for (int i=0;i<4;++i)
		{
			double ex=0;
			double ex2=0;
			double totalx=0;
			double dx=0;
			for (int j=0;j<ANGLE_LEVEL_CNT;++j)
			{
				int tv=get_block_value_int_ver(img_angle_level_block[j],r-rga[i],c-rga[i],r+rga[i],c+rga[i]);
				totalx+=tv;
			}
			if (totalx!=0)
			{
				for (int j=0;j<ANGLE_LEVEL_CNT;++j)
				{
					int tv=get_block_value_int_ver(img_angle_level_block[j],r-rga[i],c-rga[i],r+rga[i],c+rga[i]);
					ex+=tv/totalx;
					ex2+=tv*tv/totalx/totalx;
				}
				ex/=ANGLE_LEVEL_CNT;
				ex2/=ANGLE_LEVEL_CNT;
				dx=max(ex2-ex*ex,0.0);
			}
			ret.push_back((int)(dx*10000));
			if (!asn_wts) forest_wts.push_back(1);
		}
		*/
		//#define RARE_COLOR_KIND 10
		//int img_rare_level[RARE_COLOR_LEVEL][IMG_HEIGHT][IMG_WIDTH];
		//int img_rare_level_block[RARE_COLOR_LEVEL][IMG_BLOCK_HEIGHT][IMG_BLOCK_WIDTH];
		//int tmp_color_group[IMG_HEIGHT][IMG_WIDTH];
		/*
		int rgb[2]={20,40};
		for (int i=0;i<2;++i)
		{
			for (int j=0;j<RARE_COLOR_LEVEL;++j)
			{
				int rc=get_block_value_int_ver(img_rare_level_block[j],r-rgb[i],c-rgb[i],r+rgb[i],c+rgb[i]);
				ret.push_back(rc);
				if (!asn_wts) forest_wts.push_back(1);
			}
		}
		*/
		
		int rg[4]={0,15,30,40};//15,30,40
		for (int i=0;i<3;++i)
		{
			int out_tx=get_block_value_int_ver(img_grey_block[0],r-rg[i+1],c-rg[i+1],r+rg[i+1],c+rg[i+1]);
			int out_tx2=get_block_value(img_grey2_block[0],r-rg[i+1],c-rg[i+1],r+rg[i+1],c+rg[i+1]);
			int points=rg[i+1]*2+1;
			int inn_tx=0;
			int inn_tx2=0;
			if (i!=0)
			{
				inn_tx=get_block_value_int_ver(img_grey_block[0],r-rg[i],c-rg[i],r+rg[i],c+rg[i]);
				inn_tx2=get_block_value(img_grey2_block[0],r-rg[i],c-rg[i],r+rg[i],c+rg[i]);
				points-=(rg[i]*2+1);
			}
			double ex=(out_tx-inn_tx)*1.0/points;
			double ex2=(out_tx2-inn_tx2)*1.0/points;
			int dx=(int)((ex2-ex*ex)*100);
			ret.push_back(dx);
			if (!asn_wts) forest_wts.push_back(1);
						
			for (int j=0;j<GREY_LEVEL_CNT;++j)
			{
				int outer=get_block_value_int_ver(img_grey_level_block[j],r-rg[i+1],c-rg[i+1],r+rg[i+1],c+rg[i+1]);
				int inner=0;
				if (i!=0)
					inner=get_block_value_int_ver(img_grey_level_block[j],r-rg[i],c-rg[i],r+rg[i],c+rg[i]);
				
				/*
				int outer=0;
				int inner=0;
				for (int k=0;k<=j;++k)
				{
					outer+=get_block_value_int_ver(img_grey_level_block[0][k],r-rg[i+1],c-rg[i+1],r+rg[i+1],c+rg[i+1]);
					if (i!=0)
						inner+=get_block_value_int_ver(img_grey_level_block[0][k],r-rg[i],c-rg[i],r+rg[i],c+rg[i]);
				}
				*/
				ret.push_back(outer-inner);
				if (!asn_wts) forest_wts.push_back(1);
				
				
				/*
				outer=get_block_value_int_ver(img_red_level_block[0][j],r-rg[i+1],c-rg[i+1],r+rg[i+1],c+rg[i+1]);
				inner=0;
				if (i!=0)
					inner=get_block_value_int_ver(img_red_level_block[0][j],r-rg[i],c-rg[i],r+rg[i],c+rg[i]);
				ret.push_back(outer-inner);
				if (!asn_wts) forest_wts.push_back(1);
				
				outer=get_block_value_int_ver(img_green_level_block[0][j],r-rg[i+1],c-rg[i+1],r+rg[i+1],c+rg[i+1]);
				inner=0;
				if (i!=0)
					inner=get_block_value_int_ver(img_green_level_block[0][j],r-rg[i],c-rg[i],r+rg[i],c+rg[i]);
				ret.push_back(outer-inner);
				if (!asn_wts) forest_wts.push_back(1);
				
				outer=get_block_value_int_ver(img_blue_level_block[0][j],r-rg[i+1],c-rg[i+1],r+rg[i+1],c+rg[i+1]);
				inner=0;
				if (i!=0)
					inner=get_block_value_int_ver(img_blue_level_block[0][j],r-rg[i],c-rg[i],r+rg[i],c+rg[i]);
				ret.push_back(outer-inner);
				if (!asn_wts) forest_wts.push_back(1);
				*/
			}
			for (int j=0;j<GRADIENT_KIND_CNT;++j)
			{
				int outer=get_block_value_int_ver(img_gradient_level_block[j],r-rg[i+1],c-rg[i+1],r+rg[i+1],c+rg[i+1]);
				int inner=0;
				if (i!=0)
					inner=get_block_value_int_ver(img_gradient_level_block[j],r-rg[i],c-rg[i],r+rg[i],c+rg[i]);
				ret.push_back(outer-inner);
				if (!asn_wts) forest_wts.push_back(1);
			}
		}
		
		int rg2[7]={-35,-25,-15,0,15,25,36};//{-40,-30,-20,-10,0,10,20,30,41};//{-35,-20,0,20,36}//{-35,-20,0,20,36};
		for (int i=0;i<6;++i)
		{
			for (int j=0;j<GRADIENT_KIND_CNT;++j)
			{
				int v=get_block_value_int_ver(img_gradient_level_block[j],r+rg2[i],c+rg2[0],r+rg2[i+1]-1,c+rg2[6]-1);
				ret.push_back(v);
				if (!asn_wts) forest_wts.push_back(1);
			}
			for (int j=0;j<GRADIENT_KIND_CNT;++j)
			{
				int v=get_block_value_int_ver(img_gradient_level_block[j],r+rg2[0],c+rg2[i],r+rg2[6]-1,c+rg2[i+1]-1);
				ret.push_back(v);
				if (!asn_wts) forest_wts.push_back(1);
			}
			
			/*
			double points=(rg2[i+1]-rg2[i])*(rg2[6]-rg2[0]);
			int tx=get_block_value_int_ver(img_grey_block[0],r+rg2[i],c+rg2[0],r+rg2[i+1]-1,c+rg2[6]-1);
			int tx2=get_block_value(img_grey2_block[0],r+rg2[i],c+rg2[0],r+rg2[i+1]-1,c+rg2[6]-1);
			double ex=tx/points;
			double ex2=tx2/points;
			int dx=(int)((ex2-ex*ex)*100);
			ret.push_back(dx);
			if (!asn_wts) forest_wts.push_back(1);
			
			tx=get_block_value_int_ver(img_grey_block[0],r+rg2[0],c+rg2[i],r+rg2[6]-1,c+rg2[i+1]-1);
			tx2=get_block_value(img_grey2_block[0],r+rg2[0],c+rg2[i],r+rg2[6]-1,c+rg2[i+1]-1);
			ex=tx/points;
			ex2=tx2/points;
			dx=(int)((ex2-ex*ex)*100);
			ret.push_back(dx);
			if (!asn_wts) forest_wts.push_back(1);
			*/
		}
		
		int rg3=35;
		int gpid=-1;
		int gp_max_cnt=-1;
		for (int i=r-rg3;i<=r+rg3;++i)
		{
			for (int j=c-rg3;j<=c+rg3;++j)
			{
				int tid=img_gradient_indexes[i][j];
				if (tid>=0 && gra_groups[tid].cnt>gp_max_cnt)
				{
					gp_max_cnt=gra_groups[tid].cnt;
					gpid=tid;
				}
			}
		}
		int maxdif=-1;
		int mindif=-1;
		//int dif_rto=-1;
		if (gpid>=0)
		{
			maxdif=gra_groups[gpid].maxx-gra_groups[gpid].minx+1;
			mindif=gra_groups[gpid].maxy-gra_groups[gpid].miny+1;
			if (maxdif<mindif)
			{
				swap(maxdif,mindif);
			}
			/*
			if (maxdif>8)
			{
				dif_rto=maxdif*1000/mindif;
			}
			*/
		}
		ret.push_back(gp_max_cnt);
		if (!asn_wts) forest_wts.push_back(5);
		ret.push_back(maxdif);
		if (!asn_wts) forest_wts.push_back(5);
		//ret.push_back(mindif);
		//if (!asn_wts) forest_wts.push_back(1);
		
		/*
		for (int i=0;i<4;++i)
		{
			ret.push_back(img_tot_color[0][i]);
			if (!asn_wts) forest_wts.push_back(3);
		}
		*/
				
		/*
		int offr8[8]={-1,-1,-1,0,0,1,1,1};
		int offc8[8]={-1,0,1,-1,1,-1,0,1};
		struct GraGroup
		{
			int cnt;
			int maxx,minx,maxy,miny;
		};
		int now_group_index=0;
		vector<GraGroup> gra_groups;
		int img_gradient_indexes[IMG_HEIGHT][IMG_WIDTH];
		*/
		asn_wts=true;

		return ret;
	}
	void add_train_records()
	{
		int space=(int)(600/(sqrt(each_image_collect_cnt/2.0)));
		int edge_left=40;
		for (int r=edge_left;r<IMG_HEIGHT-edge_left;r+=space)
		{
			for (int c=edge_left;c<IMG_WIDTH-edge_left;c+=space)
			{
				atr_idx++;
				vector<int> prop=get_prop(r,c);
				int val=0;
				for (int i=0;i<antrec.size();++i)
				{
					int difx=antrec[i].x-c;
					int dify=antrec[i].y-r;
					int d2=difx*difx+dify*dify;
					if (antrec[i].ord<1000000&&last_mark_ord[antrec[i].ord]!=atr_idx&&d2<=1600)
					{
						last_mark_ord[antrec[i].ord]=atr_idx;
						val++;
					}
				}
				double fct=val*1.0/each_image_chk_cnt;
				fct=min(1.0,fct);
				//if(fct!=0) cerr<<r<<" "<<c<<" "<<fct<<endl;
				fact_pool.push_back(fct);
				prop_pool.push_back(prop);
			}
		}		
	}
	string receiveImage(string imageId,vector<int>& image)
	{	
		
		if (left_aloc<each_image_chk_cnt)
			return "END";
		if (each_image_chk_cnt>annotators[to_chk_img_idx].size())
		{
			to_chk_img_ord++;
			if (to_chk_img_ord>=TRAIN_IMG_CNT)
				return "END";
			to_chk_img_idx=chk_train_img_order[to_chk_img_ord].second;
			return train_img_name[to_chk_img_idx];
		}
		generate_image_info(image);
		//cerr<<imageId<<endl;
		/*
		long long chks_grey_level[8]={0,0,0,0,0,0,0,0};
		long long chks_gra_level[6]={0,0,0,0,0,0};
		for (int i=0;i<IMG_BLOCK_HEIGHT;++i)
		{
			for (int j=0;j<IMG_BLOCK_WIDTH;++j)
			{
				for (int k=0;k<8;++k)
					chks_grey_level[k]+=img_grey_level_block[1][k][i][j]*(i+j);
				for (int k=0;k<6;++k)
					chks_gra_level[k]+=img_gradient_level_block[1][k][i][j]*(i+j);
			}
		}
		*/
		//for (int i=0;i<8;++i) cerr<<chks_grey_level[i]<<" ";
		//for (int i=0;i<6;++i) cerr<<chks_gra_level[i]<<" ";
		//cerr<<endl;
				
		/*
		for(int i=0;i<GREY_LEVEL_CNT;++i)
		{
			cerr<<img_grey_level_block[0][i][IMG_HEIGHT][IMG_WIDTH]<<" ";
		}
		cerr<<endl;
		*/
		to_check_ant_idx=0;
		antrec.clear();
		antrec.reserve(each_image_chk_cnt*5);
		//until_now+=tm.read_timer();
		return imageId+","+convert_int(annotators[to_chk_img_idx][to_check_ant_idx].ord);
	}
	string receiveAnnotations(string imageId,int annotator,vector<string>& annotations)
	{
		
		for(int i=0;i<annotations.size();++i)
		{
			string str=annotations[i];
			for (int j=0;j<str.size();++j)
				if (str[j]==',')
					str[j]=' ';
			stringstream sa(str);
			AntRec ar;
			string tstr;
			sa>>tstr>>ar.x>>ar.y;
			if(tstr=="structure")
			{
				ar.ord=annotator;
				antrec.push_back(ar);
				//cmd_rec<<"fd ";
			}
		}
		--left_aloc;
		++to_check_ant_idx;
		if(to_check_ant_idx<each_image_chk_cnt)
		{
			//until_now+=tm.read_timer();
			return imageId+","+convert_int(annotators[to_chk_img_idx][to_check_ant_idx].ord);
		}
		else
		{
			add_train_records();
			++to_chk_img_ord;
			if (to_chk_img_ord<TRAIN_IMG_CNT)
			{
				to_chk_img_idx=chk_train_img_order[to_chk_img_ord].second;
				//until_now+=tm.read_timer();
				return train_img_name[to_chk_img_idx];
			}
		}
		//until_now+=tm.read_timer();
		return "END";	
	}
	vector<string> labelImage(vector<int>& image)
	{
		
		if (!build_label)
		{
			cerr<<"Before build trees seconds:"<<until_now<<endl;
			build_trees();
			//cerr<<"After build trees seconds:"<<until_now+tm.read_timer()<<endl;
			build_label=true;
		}
		generate_image_info(image);
		
		#ifdef TEST_IMG
			string fn=".\\example1_tmp\\res__"+convert_int(label_index)+"a.bmp";
			vector<int> timg;
			for (int i=0;i<IMG_HEIGHT;++i)
			{
				for (int j=0;j<IMG_WIDTH;++j)
				{
					/*
					timg.push_back((int)(nowr[tmp_color_group[i][j]]+0.5));
					timg.push_back((int)(nowg[tmp_color_group[i][j]]+0.5));
					timg.push_back((int)(nowb[tmp_color_group[i][j]]+0.5));
					*/
					
					int gry=img_grey[0][i][j];
					timg.push_back(gry);
					timg.push_back(gry);
					timg.push_back(gry);
					
					
					/*
					int tid=img_gradient_indexes[i][j];
					if (tid<0)
					{
						timg.push_back(0);
						timg.push_back(0);
						timg.push_back(0);
					}
					else
					{
						timg.push_back(r_asn[tid]);
						timg.push_back(g_asn[tid]);
						timg.push_back(b_asn[tid]);
					}
					*/
				}
			}
			gen_img(image[1],image[0],timg,fn.c_str());
		#endif
		
		vector<string> ret;
		/*
		int ii=0;
		for (int r=38;r<600;r+=69)
		{
			int cfrom=40;
			int kto=15;
			if (ii%2==1)
			{
				cfrom=81;
				kto=14;
			}
			for (int k=0;k<kto;++k)
			{
				int c=cfrom+82*k;
				vector<int> prop=get_prop(r,c);
				double pv=0;
				for (int i=0;i<TREES_CNT;++i)
					pv+=trees[i].make_decision(prop);
				pv/=TREES_CNT;
				if (pv>eps)
					ret.push_back(convert_int(c)+","+convert_int(r)+","+convert_double(pv));
			}
			ii++;
		}
		*/
		
		vector<pair<double,pair<int,int> > > cand;
		vector<pair<double,pair<int,int> > > picked;
		for (int r=40;r+40<IMG_HEIGHT;r+=6)
		{
			for (int c=40;c+40<IMG_WIDTH;c+=6)
			{
				vector<int> prop=get_prop(r,c);
				double pv=0;
				for (int i=0;i<TREES_CNT;++i)
					pv+=trees[i].make_decision(prop);
				pv/=TREES_CNT;
				if (pv>eps)
					cand.push_back(make_pair(-pv,make_pair(r,c)));
			}
		}
		sort(cand.begin(),cand.end());
		int give_con=0;
		for (int i=0;i<cand.size();++i)
		{
			bool overlap=false;
			for (int j=0;j<picked.size();++j)
			{
				int dify=cand[i].second.first-picked[j].second.first;
				int difx=cand[i].second.second-picked[j].second.second;
				if (dify*dify+difx*difx<=80*80)
				{
					overlap=true;
					break;
				}
			}
			if (!overlap)
			{
				picked.push_back(cand[i]);
				int r=cand[i].second.first;
				int c=cand[i].second.second;
				double pv=-cand[i].first;
				
				if (give_con<1)
					pv*=2.7;
				else if (give_con<2)
					pv*=2.1;
				else if (give_con<3)
					pv*=1.8;
				else if (give_con<5)
					pv*=1.7;
				else if (give_con<8)
					pv*=1.4;
				else if (give_con<12)
					pv*=1.1;
				else if (give_con<16)
					pv*=1.0;
				else if (give_con<20)
					pv*=0.9;
				else if (give_con<25)
					pv*=0.85;
				else if (give_con<30)
					pv*=0.8;
				else if (give_con<40)
					pv*=0.75;
				else if (give_con<60)
					pv*=0.7;
				else if (give_con<80)
					pv*=0.6;
				else if (give_con<100)
					pv*=0.5;
				else 
					pv*=0.3;
				
				
				pv=min(pv,1.0);
				ret.push_back(convert_int(c)+","+convert_int(r)+","+convert_double(pv));
				give_con++;
				//if(give_con==5) break;//
			}
		}
		
		label_index++;
		if (label_index==100)
		{
			//cerr<<"Tot time:"<<until_now+tm.read_timer()<<endl;
			#ifdef LOCAL
			for (int i=0;i<500;++i)
				if(maybe_cnt[i]!=0)
					cmd_rec<<"prop "<<i<<" - "<<select_cnt[i]<<"/"<<maybe_cnt[i]<<"="<<select_cnt[i]*1.0/maybe_cnt[i]*100<<endl;
			#endif
		}
		//until_now+=tm.read_timer();
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