#include <algorithm>
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <cmath>

#define AERRO 0.10		// ~ 5°
#define LTOL 100.0		// ~ 3p
#define EPS 0.00000001	// 1e-07
#define PCROSS(a1, a2, b1, b2) ((a1*b2) - (a2*b1))
#define DISTANCE(x1, y1, x2, y2) ((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2))

using namespace std;

const double PI = 2*acos(0.0);
const double PIv2 = 2*PI;
const double PId2 = PI/2;
const double PIv3d2 = (3*PI)/2;

bool between(double ref, double tst);
void rotate(double *xx, double *yy, double x, double y, double alpha);
double diff_degree(double a, double b);
double degree_vectors(double x1, double y1, double x2, double y2);
pair<double, double> build_vector(double dir);
pair<double, double> func(pair<double, double> vec_r, pair<double, double> vec_p1, pair<double, double> vec_p2);

struct Minutiae
{
	int t;//t - tipo, d - direção do vetor [0, 2*pi]
	double x, y, d;
};
struct Edge//Segmento
{
	int p1, p2;//pontos que formam o segmento
	double len_r, alphap1, alphap2;//l - tamanho do segmento, alphap1 - angulo do ponto p1, alphap2 - angulo do ponto p2
};

bool comp(Edge A, Edge B);
double match_degree(Minutiae A, Minutiae B);

struct Sample
{
	int p_size, s_size;//p_size - quantidade de pontos, s_size - quantidade de segmentos
	Minutiae *point;
	Edge *segment;
	bool p_alloc(int _n)
	{
		p_size = _n;
		point = (Minutiae*) malloc (_n*sizeof(Minutiae));
		return point != NULL;
	}
	bool s_alloc(int _n)
	{
		s_size = ((_n-1)*(_n-1) + (_n-1))>>1;
		segment = (Edge*) malloc (s_size*sizeof(Edge));
		return segment != NULL;
	}
	void p_dalloc(void)
	{
		free(point);
	}
	void s_dalloc(void)
	{
		free(segment);
	}
	void s_sort(void)
	{
		sort(segment, segment + s_size, comp);
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

				pair<double, double> ans = func(build_vector(atan2(point[j].y - point[i].y, point[j].x - point[i].x)),
													vec_p1, build_vector(point[j].d));

				segment[k].alphap1 = ans.first;
				segment[k].alphap2 = ans.second;
				k++;
			}
		}
	}
};

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
	fout = fopen(argv[3], "a+");
	if(fout == NULL)
	{
		printf("Erro ao abrir %s.\n", argv[3]);
		return 0;
	}

	int n1, n2;
	fscanf(fsp1, "%d", &n1);
	fscanf(fsp2, "%d", &n2);

	if(!A.p_alloc(n1) || !A.s_alloc(n1) || !B.p_alloc(n2) || !B.s_alloc(n2) || !M.p_alloc(n2))
	{
		puts("Erro na alocacao de memoria.");
		return 0;
	}

	A.p_read(fsp1);//Ref Sample
	B.p_read(fsp2);//Query Sample

	A.generate_segments();//PEGUE OS TIPOS TAMBÉM, PARA DEIXAR MAIS CONFIÁVEL!!!
	B.generate_segments();

	A.s_sort();
	B.s_sort();

	int l = 0, r = 0, cont = 0;
	double max_sum = 0;

	for(int i = 0;i < A.s_size;i++)
	{//[l, r)

		double ll = A.segment[i].len_r - LTOL, rr = A.segment[i].len_r + LTOL;

		while(B.segment[l].len_r < ll && !(l > B.s_size))
			l++;
		while(B.segment[r].len_r < rr && !(r > B.s_size))
			r++;

		int id_ap1 = A.segment[i].p1, id_ap2 = A.segment[i].p2;

		for(int j = l;j < r;j++)
		{
			int id_bp1 = B.segment[j].p1, id_bp2 = B.segment[j].p2;
			
			if(between(A.segment[i].alphap1, B.segment[j].alphap1) && between(A.segment[i].alphap2, B.segment[j].alphap2))
				if(A.point[id_ap1].t == B.point[id_bp1].t && A.point[id_ap2].t == B.point[id_bp2].t)//Os tipos tem que ser iguais
				{
					cont++;
					double theta = atan2(A.point[id_ap2].y - A.point[id_ap1].y, A.point[id_ap2].x - A.point[id_ap1].x)
									- atan2(B.point[id_bp2].y - B.point[id_bp1].y, B.point[id_bp2].x - B.point[id_bp1].x);

					for(int k = 0;k < B.p_size;k++)
					{
						rotate(&M.point[k].x, &M.point[k].y, B.point[k].x, B.point[k].y, theta);

						M.point[k].d = B.point[k].d + theta;
						if(M.point[k].d < 0)
							M.point[k].d += PIv2;
						if(M.point[k].d > PIv2)
							M.point[k].d -= PIv2;

						M.point[k].t = B.point[k].t;
					}

					double hh = ((A.point[id_ap1].x - M.point[id_bp1].x) + (A.point[id_ap2].x - M.point[id_bp2].x))/2;
					double kk = ((A.point[id_ap1].y - M.point[id_bp1].y) + (A.point[id_ap2].y - M.point[id_bp2].y))/2;

					bool flag = true;
					for(int k = 0;k < B.p_size && flag;k++)
					{
						M.point[k].x += hh;
						M.point[k].y += kk;
						if(M.point[k].x < 0.0 || M.point[k].y < 0.0 || M.point[k].x > 500 || M.point[k].y > 500)
							flag = false;
					}
					if(!flag)
						continue;

					//printf("%d\n", M.p_size);
					//M.p_print();

					double sum = 0;
					for(int k = 0;k < A.p_size;k++)
						for(int pp = 0;pp < B.p_size;pp++)
							sum += match_degree(A.point[k], M.point[pp]);

					max_sum = max(max_sum, sum);
				}
		}
	}

	double anss = (max_sum/max(n1, n2))*100.00;
	if(anss > 100)
		anss = 100;
    int ass = (int)round(anss);
    if(ass == 0)
        fprintf(fout, "%s;%s;FAIL;0\n", argv[1], argv[2]);
    else
        fprintf(fout, "%s;%s;OK;%d\n", argv[1], argv[2], ass);
	//printf("%d %lf\n", cont, (max_sum/max(n1, n2))*100.00);

	A.p_dalloc();
	A.s_dalloc();
	B.p_dalloc();
	B.s_dalloc();
	M.p_dalloc();

	fclose(fsp1);
	fclose(fsp2);
	fclose(fout);
	return 0;
}
double degree_vectors(double x1, double y1, double x2, double y2)//Calcula o angulo entre dois vetores
{
	return acos((x1*x2 + y1*y2)/(sqrt(DISTANCE(0, 0, x1, y1))*sqrt(DISTANCE(0, 0, x2, y2))));
}
pair<double, double> build_vector(double dir)
{
	double x = 1, y = tan(dir);
	if(fabs(y - 0.0) < EPS)
	{
		if(dir < PId2)
			return make_pair(x, 0);
		else if(dir < PI)
			return make_pair(-x, 0);
		else if(dir < PIv3d2)
			return make_pair(-x, 0);
		else
			return make_pair(x, 0);
	}
	else
	{
		if(dir < PId2)
			return make_pair(x, y);
		else if(dir < PI)
			return make_pair(-x, -y);
		else if(dir < PIv3d2)
			return make_pair(-x, -y);
		else
			return make_pair(x, y);
	}
}
pair<double, double> func(pair<double, double> vec_r, pair<double, double> vec_p1, pair<double, double> vec_p2)
{	//p1 sempre abaixo de p2
	if(vec_r.second < 0)
		return make_pair(degree_vectors(-1*vec_r.first, -1*vec_r.second, vec_p1.first, vec_p1.second),
							degree_vectors(vec_r.first, vec_r.second, vec_p2.first, vec_p2.second));
	else
		return make_pair(degree_vectors(vec_r.first, vec_r.second, vec_p1.first, vec_p1.second),
							degree_vectors(-1*vec_r.first, -1*vec_r.second, vec_p2.first, vec_p2.second));
}
bool comp(Edge A, Edge B)
{
	if(A.len_r == B.len_r)
	{
		if(A.alphap1 == B.alphap1)
		{
			if(A.alphap2 == B.alphap2)
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
		if(tst < (PIv2 + ll) && tst > rr)
			return false;
		else
			return true;
	}
	else if(rr > PIv2)
	{
		if(tst < ll && tst > ((rr - PIv2)))
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
void rotate(double *xx, double *yy, double x, double y, double alpha)
{
	*(xx) = (x*cos(alpha)) - (y*sin(alpha));
	*(yy) = (x*sin(alpha)) + (y*cos(alpha));
}
double match_degree(Minutiae A, Minutiae B)//tentar fazer um match linear
{
	double score = 0, dt = DISTANCE(A.x, A.y, B.x, B.y), diff = diff_degree(A.d, B.d);

	if(dt < 400.0)
		score += (400.0 - dt)/400.0;
	else
		return 0;

	if(diff < 0.20)
		score += (0.20 - diff)/0.20;
	else if(diff > (PI - 0.20))
		score += (0.20 - (PI - diff))/0.20;
	else
		return 0;

	if(A.t == B.t)
		score += score/2;
	else
		score += score/4;

	return score/3;
}
double diff_degree(double a, double b)//diferença entre dois angulos
{
	double d = fabs(a - b);
	return min(d, PIv2 - d);
}
