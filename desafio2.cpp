#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <cstdio>
#include <cmath>

#define PCROSS(a1, a2, b1, b2) ((a1*b2) - (a2*b1))
#define DISTANCE(x1, y1, x2, y2) ((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2))

using namespace std;

const double PI = 2*acos(0.0);
const double PIv2 = 2*PI;
const double PId2 = PI/2;
const double PIv3d2 = (3*PI)/2;

struct Minutiae
{
	int t;//t - tipo, d - direção do vetor [0, 2*pi]
	double x, y, d;
};
struct Sample
{
	int p_size;//p_size - quantidade de pontos, s_size - quantidade de segmentos
	Minutiae *point;
	bool p_alloc(int _n)
	{
		p_size = _n;
		point = (Minutiae*) malloc (_n*sizeof(Minutiae));
		return point != NULL;
	}
	void p_dalloc(void)
	{
		free(point);
	}
	void p_read(FILE *file)
	{
		for(int i = 0;i < p_size;i++)
			fscanf(file, "%lf %lf %lf %d", &point[i].x, &point[i].y, &point[i].d, &point[i].t);
	}
	void p_print(void)
	{
		for(int i = 0;i < p_size;i++)
			printf("%lf %lf %lf %d\n", point[i].x, point[i].y, point[i].d, point[i].t);
	}
};

void priority_process(int pq[], int id, Sample *A);
void rotate(double *xx, double *yy, double x, double y, double alpha);
double diff_degree(double a, double b);
double match_degree(Minutiae A, Minutiae B);

int main(int argc, char **argv)
{
	FILE *fsp1, *fsp2, *fout;
	Sample A, B, M;
	
	fsp1 = fopen(argv[1], "r");
	if(fsp1 == NULL)
	{
		printf("Erro ao abrir %s.\n", argv[1]);
		return 0;
	}
	fsp2 = fopen(argv[2], "r");
	if(fsp2 == NULL)
	{
		printf("Erro ao abrir %s.\n", argv[2]);
		return 0;
	}
	fout = fopen(argv[3], "w");
	if(fout == NULL)
	{
		printf("Erro ao abrir %s.\n", argv[3]);
		return 0;
	}
	
	int n1, n2;
	fscanf(fsp1, "%d", &n1);
	fscanf(fsp2, "%d", &n2);
	
	if(!A.p_alloc(n1) || !B.p_alloc(n2) || !M.p_alloc(n2))
	{
		puts("Erro na alocacao de memoria.");
		return 0;
	}
	
	A.p_read(fsp1);//Ref Sample
	B.p_read(fsp2);//Query Sample
	
	vector<int>hash[35][35], bestinM[n1];
	int pq[n1];
	int combi[9] = {0,  0, 0, -1, 1, 1,  1, -1, -1};
	int combj[9] = {0, -1, 1,  0, 0, 1, -1,  1, -1};
	
	for(int i = 0;i < n1;i++)
	{
		int j = A.point[i].x/15.0;
		int k = A.point[i].y/15.0;
		hash[j][k].push_back(i);
	}
	
	double max_sum = 0;
	
	for(int i = 0;i < n1;i++)
	{
		priority_process(pq, i, &A);
		for(int j = 0;j < n2;j++)
		{
			double theta = A.point[i].d - B.point[j].d;
			for(int k = 0;k < n2;k++)
			{
				rotate(&M.point[k].x, &M.point[k].y, B.point[k].x, B.point[k].y, theta);
					
				M.point[k].d = B.point[k].d + theta;
				if(M.point[k].d < 0)
					M.point[k].d += PIv2;
				if(M.point[k].d > PIv2)
					M.point[k].d -= PIv2;
					
				M.point[k].t = B.point[k].t;
			}
			
			double hh = A.point[i].x - M.point[j].x;
			double kk = A.point[i].y - M.point[j].y;
			for(int k = 0;k < n2;k++)
			{
				M.point[k].x += hh;
				M.point[k].y += kk;
			}
			
			for(int k = 0;k < n2;k++)
			{	//os melhores para cada ponto de M
				int p = M.point[k].x/15.0;
				int q = M.point[k].y/15.0;
				
				for(int ii = 0;ii < 9;ii++)
				{
					int jj = p + combi[ii];
					int kk = q + combj[ii];
					if(kk >= 0 && jj >= 0 && jj < 35 && kk < 35)//EXCLUI O MATCHING DE PONTOS NEGATIVOS OU ACIMA DE 500
						for(int idx : hash[jj][kk])
							if(DISTANCE(A.point[idx].x, A.point[idx].y, M.point[k].x, M.point[k].y) < 225.0)//testar aumentar e diminuir
								bestinM[idx].push_back(k);
				}
			}
			
			int cont = 0;
			for(int k = 0;k < n1;k++)
				if(!bestinM[k].empty())
					cont++;
			
			vector<bool>used(n2, false);
			double sum = 0;
			for(int k = 0;k < n1;k++)
			{
				int idd = -1;
				double mx = -1;
				for(int idx : bestinM[pq[k]])//best in M to p[k]
					if(!used[idx])
					{
						double score = match_degree(A.point[pq[k]], M.point[idx]);
						if(score > mx)
						{
							mx = score;
							idd = idx;
						}
					}
				
				if(idd != -1)
				{
					sum += mx;
					used[idd] = true;
				}
			}
				
			max_sum = max(max_sum, sum/cont);
			
			for(int k = 0;k < n1;k++)
					bestinM[k].clear();
		}
	}
		
	printf("%lf\n", max_sum*100.00);
	
	A.p_dalloc();
	B.p_dalloc();
	M.p_dalloc();
	
	fclose(fsp1);
	fclose(fsp2);
	fclose(fout);
	return 0;
}
void rotate(double *xx, double *yy, double x, double y, double alpha)
{
	*(xx) = (x*cos(alpha)) - (y*sin(alpha));
	*(yy) = (x*sin(alpha)) + (y*cos(alpha));
}
double match_degree(Minutiae A, Minutiae B)//tentar fazer um match linear
{
	double score = 0, dt = DISTANCE(A.x, A.y, B.x, B.y), diff = diff_degree(A.d, B.d);
	
	if(dt < 625)
		score += (625 - dt)/625.0;
	else
		return 0;
		
	if(diff < 0.25)
		score += (0.25 - diff)/0.25;
	else if(diff > (PI - 0.1))
		score += (0.1 - (PI - diff))/0.1;
	else
		return 0;
		
	if(A.t == B.t)
		score += score/2;
	else
		score += score/8;
		
	return score/3;
}
double diff_degree(double a, double b)//diferença entre dois angulos
{
	double d = fabs(a - b);
	return min(d, PIv2 - d);
}
void priority_process(int pq[], int id, Sample *A)
{
	vector<pair<double, int> >val(A->p_size);
	double x = A->point[id].x, y = A->point[id].y;
	
	for(int i = 0;i < A->p_size;i++)
		val[i] = make_pair(DISTANCE(x, y, A->point[i].x, A->point[i].y), i);
								
	sort(val.begin(), val.end());
	
	for(int i = 0;i < val.size();i++)
		pq[i] = val[i].second;
}
