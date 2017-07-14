#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#define AERRO 0.1		// ~ 5°
#define LTOL 1.5		// ~ sqrt(2)
#define EPS 0.00000001	// 1e-07
#define PCROSS(a1, a2, b1, b2) ((a1*b2) - (a2*b1))
#define DISTANCE(x1, y1, x2, y2) (sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2)))
using namespace std;
const double PI = 2*acos(0.0);

bool between(double ref, double tst);
double degree_vectors(double x1, double y1, double x2, double y2);
pair<double, double> build_vector(double dir);
pair<double, double> func(bool op, pair<double, double> vec_r, pair<double, double> vec_p1, pair<double, double> vec_p2);

struct Minutiae
{
	int x, y, t;//t - tipo, d - direção do vetor [0, 2*pi]
	double d;
};
struct Edge//Segmento
{
	int p1, p2;//pontos que formam o segmento
	double len_r, alphap1, alphap2;//l - tamanho do segmento, cf_r - coeficiente do segmento, alphap1 - angulo do ponto p1 com o segmento, alphap2 - angulo do ponto p2 com o segmento
};

bool comp(Edge A, Edge B);

struct Sample
{
	int p_size, s_size;//p_size - quantidade de pontos, s_size - quantidade de segmentos
	Minutiae *point;
	Edge *segment;
	bool p_alloc(int _n)
	{
		p_size = _n;
		s_size = ((_n-1)*(_n-1) + (_n-1))>>1;
		point = (Minutiae*) malloc (_n*sizeof(Minutiae));
		if(point == NULL)
			return false;
		segment = (Edge*) malloc (s_size*sizeof(Edge));
		if(segment == NULL)
			return false;
		return true;
	}
	void p_dalloc(void)
	{
		free(point);
		free(segment);
	}
	void p_read(FILE *file)
	{
		for(int i = 0;i < p_size;i++)
			fscanf(file, "%d %d %lf %d", &point[i].x, &point[i].y, &point[i].d, &point[i].t);
	}
	void p_print(void)
	{
		for(int i = 0;i < p_size;i++)
			printf("%d %d %lf %d -> %lf\n", point[i].x, point[i].y, point[i].d, point[i].t, tan(point[i].d));
	}
	void s_print(void)
	{
		for(int i = 0;i < s_size;i++)
			printf("(%d, %d) -> %lf %lf %lf\n", segment[i].p1, segment[i].p2, segment[i].len_r, segment[i].alphap1, segment[i].alphap2);
	}
	void generate_segments(void)
	{
		for(int i = 0, k = 0;i < p_size;i++)
		{
			pair<double, double>vec_p1 = build_vector(point[i].d);
			for(int j = i+1;j < p_size;j++)
			{
				segment[k].p1 = i;
				segment[k].p2 = j;
				segment[k].len_r = DISTANCE(point[i].x, point[i].y, point[j].x, point[j].y);
				
				pair<double, double> ans;
				if(point[i].x != point[j].x)
				{
					if(point[i].y < point[j].y || (point[i].y == point[j].y && point[i].x < point[j].x))
						ans = func(true, build_vector(atan2((double)(point[j].y - point[i].y), (double)(point[j].x - point[i].x))), vec_p1, build_vector(point[j].d));
					else
						ans = func(false, build_vector(atan2((double)(point[j].y - point[i].y), (double)(point[j].x - point[i].x))), vec_p1, build_vector(point[j].d));
				}
				else
				{
					if(point[i].y < point[j].y && (point[i].y == point[j].y && point[i].x < point[j].x))
						ans = func(true, make_pair(0, 1), vec_p1, build_vector(point[j].d));
					else
						ans = func(false, make_pair(0, 1), vec_p1, build_vector(point[j].d));
				}
				segment[k].alphap1 = ans.first;
				segment[k].alphap2 = ans.second;
				
				k++;
			}
		}
	}
	void s_sort(void)
	{
		sort(segment, segment + s_size, comp);
	}
};
int main(int argc, char **argv)
{
	FILE *fsp1, *fsp2, *fout;
	Sample A, B;
	
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
	
	if(!A.p_alloc(n1) || !B.p_alloc(n2))
	{
		puts("Erro na alocacao de memoria.");
		return 0;
	}
	
	printf("%d %d\n%d %d\n", A.p_size, B.p_size, A.s_size, B.s_size);
	
	A.p_read(fsp1);
	B.p_read(fsp2);
	
	A.generate_segments();
	B.generate_segments();
	
	A.s_sort();
	B.s_sort();
	
	int l = 0, r = 0, cont = 0;
	for(int i = 0;i < A.s_size;i++)
	{
		//[l, r)
		double ll = A.segment[i].len_r - LTOL, rr = A.segment[i].len_r + LTOL;
		while(B.segment[l].len_r < ll && !(l > B.s_size))
			l++;
		while(B.segment[r].len_r < rr && !(r > B.s_size))
			r++;
		
		int id_ap1 = A.segment[i].p1, id_ap2 = A.segment[i].p2;
		
		for(int j = l;j < r;j++)
		{						//p1 com p1, p2 com p2
			int id_bp1 = B.segment[j].p1, id_bp2 = B.segment[j].p2;
			if(between(A.segment[i].alphap1, B.segment[j].alphap1) && between(A.segment[i].alphap2, B.segment[j].alphap2))
			{
				//Ap1 -> Ap2 and Bp1 -> Bp2	
				pair<int, int> vec_a = make_pair(A.point[id_ap2].x - A.point[id_ap1].x, A.point[id_ap2].y - A.point[id_ap1].y);
				pair<int, int> vec_b = make_pair(B.point[id_bp2].x - B.point[id_bp1].x, B.point[id_bp2].y - B.point[id_bp1].y);
				
				if(PCROSS(vec_a.first, vec_a.second, vec_b.first, vec_b.second) > 0)
				{
					//a -> b => Ap2 go to 
				}
				else
				{
					//b -> a
				}
				cont++;
			}
			else if(between(A.segment[i].alphap1, B.segment[j].alphap2) && between(A.segment[i].alphap2, B.segment[j].alphap1))
			{
				//Ap1 -> Ap2 and //Bp2 -> Bp1
				pair<int, int> vec_a = make_pair(A.point[id_ap2].x - A.point[id_ap1].x, A.point[id_ap2].y - A.point[id_ap1].y);
				pair<int, int> vec_b = make_pair(B.point[id_bp1].x - B.point[id_bp2].x, B.point[id_bp1].y - B.point[id_bp2].y);
				
				if(PCROSS(vec_a.first, vec_a.second, vec_b.first, vec_b.second) > 0)
				{
					
				}
				else
				{
					
				}
				cont++;
			}
		}
	}
	
	printf("%d\n", cont);
	
	A.p_dalloc();
	B.p_dalloc();
	
	fclose(fsp1);
	fclose(fsp2);
	fclose(fout);
	return 0;
}
double degree_vectors(double x1, double y1, double x2, double y2)//Calcula o angulo entre dois vetores
{
	//printf("::\t%lf %lf %lf %lf | %lf / (%lf * %lf) | %lf = %lf\n", x1, y1, x2, y2, x1*x2 + y1*y2, DISTANCE(0, 0, x1, y1), DISTANCE(0, 0, x2, y2), DISTANCE(0, 0, x1, y1)*DISTANCE(0, 0, x2, y2), acos((x1*x2 + y1*y2)/(DISTANCE(0, 0, x1, y1)*DISTANCE(0, 0, x2, y2))));
	return acos((x1*x2 + y1*y2)/(DISTANCE(0, 0, x1, y1)*DISTANCE(0, 0, x2, y2)));
}
pair<double, double> build_vector(double dir)
{
	double x = 1, y = tan(dir);
	if(fabs(y - 0.0) < EPS)
	{
		if(dir < PI/2)
			return make_pair(x, 0);
		else if(dir < PI)
			return make_pair(-x, 0);
		else if(dir < (3*PI)/2)
			return make_pair(-x, 0);
		else
			return make_pair(x, 0);	
	}
	else
	{
		if(dir < PI/2)
			return make_pair(x, y);
		else if(dir < PI)
			return make_pair(-x, -y);
		else if(dir < (3*PI)/2)
			return make_pair(-x, -y);
		else
			return make_pair(x, y);	
	}
}
pair<double, double> func(bool op, pair<double, double> vec_r, pair<double, double> vec_p1, pair<double, double> vec_p2)
{//rever caso dy = 0;
	//printf("%d %lf %lf %lf %lf %lf %lf\n", op, vec_r.first, vec_r.second, vec_p1.first, vec_p1.second, vec_p2.first, vec_p2.second);
	//cout<<op<<endl;
	if(op)
	{
		if(vec_r.second < 0)
			return make_pair(degree_vectors(-1*vec_r.first, -1*vec_r.second, vec_p1.first, vec_p1.second), 
								degree_vectors(vec_r.first, vec_r.second, vec_p2.first, vec_p2.second));
		else
			return make_pair(degree_vectors(vec_r.first, vec_r.second, vec_p1.first, vec_p1.second), 
								degree_vectors(-1*vec_r.first, -1*vec_r.second, vec_p2.first, vec_p2.second));
	}
	else
	{
		if(vec_r.second < 0)
			return make_pair(degree_vectors(vec_r.first, vec_r.second, vec_p1.first, vec_p1.second), 
								degree_vectors(-1*vec_r.first, -1*vec_r.second, vec_p2.first, vec_p2.second));
		else
			return make_pair(degree_vectors(-1*vec_r.first, -1*vec_r.second, vec_p1.first, vec_p1.second), 
								degree_vectors(vec_r.first, vec_r.second, vec_p2.first, vec_p2.second));
	}
}
bool comp(Edge A, Edge B)
{
	if(fabs(A.len_r - B.len_r) < EPS)
	{
		if(fabs(A.alphap1 - B.alphap1) < EPS)
		{
			if(fabs(A.alphap2 - B.alphap2) < EPS)
			{
				if(A.p1 == B.p1)
					return A.p2 < B.p2;
				return A.p1 < B.p1;
			}
			return A.alphap2 < B.alphap2;
		}
		return A.alphap1 < B.alphap2;
	}
	return A.len_r < B.len_r;
}
bool between(double ref, double tst)
{//tst > ref - error && tst < ref + error
	double ll = ref - AERRO, rr = ref + AERRO;
	if(ll < 0)
	{
		if(tst < ((2*PI) + ll) && tst > rr)
			return false;
		else
			return true;
	}
	else if(rr > (2*PI))
	{
		if(tst < ll && tst > ((rr - (2*PI))))
			return false;
		else
			return true;
	}
	else
	{
		if(tst > ll && tst < rr)
			return true;
		else
			return false;
	}
}
