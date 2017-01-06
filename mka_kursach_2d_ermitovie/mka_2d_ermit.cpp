#define _CRT_SECURE_NO_WARNINGS

#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <set>
#include <locale>
#include "glut.h"
using namespace std;

ofstream outFile("output.txt");
bool StepLikaTelma = true,
	DrowText = false,
	RandomColors = false,
	OnlyWhiteColor = false;

struct colour // цвет точки
{
	int red;
	int green;
	int blue;
};

struct point // точка
{
	point(double newX, double newY)
	{
		x = newX;
		y = newY;
	}

	point();
	double x, y;

	friend bool operator==(const point& lhs, const point& rhs)
	{
		return lhs.x == rhs.x
			&& lhs.y == rhs.y;
	}

	friend bool operator!=(const point& lhs, const point& rhs)
	{
		return !(lhs == rhs);
	}
};

struct locateOfPoint
{
	int i, j;

	locateOfPoint(int new_i, int new_j)
	{
		i = new_i;
		j = new_j;
	}

	locateOfPoint()
	{
		i = j = 0;
	}
};

struct field
{
	double x1, x2, y1, y2, lambda, gamma;

	field(double x1_new, double x2_new, double y1_new, double y2_new, double lambda_new, double gamma_new)
	{
		x1 = x1_new;
		x2 = x2_new;
		y1 = y1_new;
		y2 = y2_new;
		lambda = lambda_new;
		gamma = gamma_new;
	}

	field()
	{
		x1 = x2 = y1 = y2 = lambda = gamma = 0;
	}
};

struct nvtr
{
	nvtr(int uz1, int uz2, int uz3, int uz4, int numberField_new)
	{
		uzel[0] = uz1;
		uzel[1] = uz2;
		uzel[2] = uz3;
		uzel[3] = uz4;
		numberField = numberField_new;
	}

	nvtr(int uz1, int uz2, int uz3, int uz4)
	{
		uzel[0] = uz1;
		uzel[1] = uz2;
		uzel[2] = uz3;
		uzel[3] = uz4;
		numberField = 0;
	}

	nvtr()
	{
		uzel[0] = uzel[1] = uzel[2] = uzel[3] = numberField = 0;
	}

	int uzel[4], numberField;
};

vector<point> xy;
vector<nvtr> KE;
vector<field> sreda;

GLint Width, Height;
int Red = 255, Green = 255, Blue = 255; // цвет по умолчанию

////Демонстрация различия согласованных и не согласованных
int nX, nY;
double hx = 0.01, hy = 0.01, mnojX = 1.3, mnojY = 1.3;
double leftX = 5, rightX = 25;
double leftY = 5, rightY = 25;

double* xNet;
double* yNet;

int countOfNodes;
int countOfGlobalBasicFunc;

//matrix
int *ig, *jg;
double *ggl, *ggu, *di, *b, *q;
int iter;

double localB[16], localMatrix[16][16];

double hForProizv = 1e-6;

void Reshape(GLint w, GLint h) // При изменении размеров окна
{
	Width = w;
	Height = h;
	glViewport(0, 0, w, h);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(0, w, 0, h, -1.0, 1.0);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}

int FindNumberOfKe(int uzel1, int uzel2)
{
	int i, j, k;
	for (i = 0; i < KE.size(); i++)
	{
		for (j = 0; j < 4; j++)
		{
			if (KE[i].uzel[j] == uzel1 || KE[i].uzel[j] == uzel2)
				for (k = j + 1; k < 4; k++)
				{
					if (KE[i].uzel[k] == uzel1 || KE[i].uzel[k] == uzel2)
						return i;
				}
		}
	}
	cout << "Ошибка. Не найден КЭ." << endl;
	system("pause");
	exit(1);
}

int indexXY(point goal)
{
	for (int i = 0; i < xy.size(); i++)
	{
		if (goal == xy[i])
			return i;
	}

	cout << "Ошибка. Не найдена точка (" << goal.x << "," << goal.y << ")" << endl;
	system("pause");
	exit(1);
}

int indexXY(double x, double y)
{
	point goal(x, y);
	return indexXY(goal);
}

int NumberNode(double x, double y)
{
	int i;
	for (i = 0; i < xy.size(); i++)
	{
		if (xy[i].x == x && xy[i].y == y)
			return i;
	}
	return -1;
}

locateOfPoint FindLocate(point sample)
{
	locateOfPoint a;
	a.i = -1;
	a.j = -1;
	int i, j;
	for (i = 0; i < nX; i++)
	{
		if (xNet[i] == sample.x)
			a.i = i;
	}
	for (j = 0; j < nY; j++)
	{
		if (yNet[j] == sample.y)
			a.j = j;
	}
	if (a.i == -1 || a.j == -1)
	{
		cout << "Ошибка FindLocate: не найдена точка" << endl;
		system("pause");
		exit(1);
	}
	return a;
}

int FindLocate(double* massiv, int razm, double x)
{
	int i;
	for (i = 0; i < razm; i++)
	{
		if (massiv[i] == x)
			return i;
	}
	return -1;
}

void drawText(const char* text, int length, int x, int y)
{
	glMatrixMode(GL_PROJECTION);
	double* matrix = new double[16];
	glGetDoublev(GL_PROJECTION_MATRIX, matrix);
	glLoadIdentity();
	glOrtho(0, 600, 0, 600, -5, 5);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glPushMatrix();
	glLoadIdentity();
	glRasterPos2i(x, y);
	for (int i = 0; i < length; i++)
	{
		glutBitmapCharacter(GLUT_BITMAP_9_BY_15, (int)text[i]);//GLUT_BITMAP_9_BY_15
	}
	glPopMatrix();
	glMatrixMode(GL_PROJECTION);
	glLoadMatrixd(matrix);
	glMatrixMode(GL_MODELVIEW);
}

void setColour(int col) // преобразование кодов цветов
{
	if (col == 7) //Каждый
	{
		Red = 255;
		Green = 0;
		Blue = 0;
	}
	if (col == 3) //Охотник
	{
		Red = 255;
		Green = 102;
		Blue = 0;
	}
	if (col == 2) //Желает
	{
		Red = 255;
		Green = 255;
		Blue = 0;
	}
	if (col == 1) //Знать
	{
		Red = 0;
		Green = 255;
		Blue = 0;
	}
	if (col == 5) //Где
	{
		Red = 0;
		Green = 255;
		Blue = 255;
	}
	if (col == 6) //Сидит
	{
		Red = 0;
		Green = 0;
		Blue = 255;
	}
	if (col == 4) //Фазан
	{
		Red = 128;
		Green = 0;
		Blue = 255;
	}
	if (col == 0) //Белый
	{
		Red = 255;
		Green = 255;
		Blue = 255;
	}
}

void Display(void) // функция вывода
{
	int i, j, k, t, colorSreda;
	glClearColor(1, 1, 1, 1); // очистка буфера
	glClear(GL_COLOR_BUFFER_BIT);
	glColor3ub(0, 0, 0);
	glShadeModel(GL_FLAT); // отключение интерполяции
	double h_x = xNet[nX - 1] - xNet[0];
	double h_y = yNet[nY - 1] - yNet[0];
	double E_x = -xNet[0];
	double E_y = -yNet[0];
	double tempY;
	if (h_x < h_y)tempY = 585 / h_y;
	else tempY = 585 / h_x;

	glLineWidth(2);
	colorSreda = 0;
	for (i = 0; i < sreda.size(); i++)
	{
		if (!RandomColors)
			colorSreda = sreda[i].gamma;
		else
		if (colorSreda > 7) colorSreda = 0;
		if (OnlyWhiteColor) colorSreda = 0;
		setColour(colorSreda);
		colorSreda++;
		glColor3ub(Red, Green, Blue);

		glBegin(GL_QUADS);
		glVertex2f((E_x + sreda[i].x1) * tempY, (E_y + sreda[i].y1) * tempY);
		glVertex2f((E_x + sreda[i].x2) * tempY, (E_y + sreda[i].y1) * tempY);
		glVertex2f((E_x + sreda[i].x2) * tempY, (E_y + sreda[i].y2) * tempY);
		glVertex2f((E_x + sreda[i].x1) * tempY, (E_y + sreda[i].y2) * tempY);
		glEnd();
		glColor3ub(0, 0, 0);
		glBegin(GL_LINE_LOOP);
		glVertex2f((E_x + sreda[i].x1) * tempY, (E_y + sreda[i].y1) * tempY);
		glVertex2f((E_x + sreda[i].x2) * tempY, (E_y + sreda[i].y1) * tempY);
		glVertex2f((E_x + sreda[i].x2) * tempY, (E_y + sreda[i].y2) * tempY);
		glVertex2f((E_x + sreda[i].x1) * tempY, (E_y + sreda[i].y2) * tempY);
		glEnd();
	}

	glLineWidth(1);
	for (j = 0; j < nY - 1; j++)
	{
		for (i = 0; i < nX - 1; i++)
		{
			glBegin(GL_LINE_LOOP);
			glVertex2f((E_x + xNet[i]) * tempY, (E_y + yNet[j]) * tempY);
			glVertex2f((E_x + xNet[i + 1]) * tempY, (E_y + yNet[j]) * tempY);
			glVertex2f((E_x + xNet[i + 1]) * tempY, (E_y + yNet[j + 1]) * tempY);
			glVertex2f((E_x + xNet[i]) * tempY, (E_y + yNet[j + 1]) * tempY);
			glEnd();
		}
	}
	glColor3ub(255, 0, 0);
	glPointSize(2);
	xy.size();
	glBegin(GL_POINTS);
	for (i = 0; i < nX; i++)
	{
		for (j = 0; j < nY; j++)
		{
			glVertex2f((E_x + xNet[i]) * tempY, (E_y + yNet[j]) * tempY);
		}
	}
	glEnd();

	glFinish();
}

void GenerateNetLikeTelma(set<double>& mas, ifstream& fileNet)
{
	double firstElement, startPosition, endPosition, position, lastPosition, stepReal;;
	int i, numberOfInterval;
	fileNet >> firstElement >> numberOfInterval;
	double* intervals = new double[numberOfInterval + 1];
	double* sizeOfSteps = new double[numberOfInterval];
	double* mnojiteli = new double[numberOfInterval];
	int* napravlenie = new int[numberOfInterval];
	intervals[0] = firstElement;
	for (i = 1; i < numberOfInterval + 1; i++)
	{
		fileNet >> intervals[i];
	}
	for (i = 0; i < numberOfInterval; i++)
	{
		fileNet >> sizeOfSteps[i];
	}
	for (i = 0; i < numberOfInterval; i++)
	{
		fileNet >> mnojiteli[i];
		if (mnojiteli[i] < 1)
		{
			cout << "Ошибка описания сетки по оси. Коэфициент разрядки должен быть больше или равен 1";
			system("pause");
			exit(1);
		}
	}
	for (i = 0; i < numberOfInterval; i++)
	{
		fileNet >> napravlenie[i];
	}
	for (i = 0; i < numberOfInterval; i++)
	{
		if (napravlenie[i] == -1)
		{
			startPosition = intervals[i + 1];
			endPosition = intervals[i];
		}
		else
		{
			startPosition = intervals[i];
			endPosition = intervals[i + 1];
		}
		stepReal = 0;
		lastPosition = startPosition;
		position = startPosition + napravlenie[i] * sizeOfSteps[i];
		while (position * napravlenie[i] < endPosition * napravlenie[i])
		{
			mas.insert(position);
			stepReal = fabs(lastPosition - position) * mnojiteli[i];
			lastPosition = position;
			position += napravlenie[i] * stepReal;
		}
		if (i != numberOfInterval - 1)
		{
			mas.insert(endPosition);
		}
		//if (fabs(lastPosition - endPosition) < sizeOfSteps[i] && lastPosition != startPosition)
		//{
		//	mas.erase(mas.find(lastPosition));
		//}
	}
}

//V - вершины
//R - регулярные узлы, 4 ребра примыкают
//Y - отсутствует вершина, 2 ребра

//W - сверху нет ребра
//A - слева нет ребра
//S - снизу нет ребра
//D - справа нет ребра
int GenerateNet()
{
	int i, numberFields;
	field tempField;
	set<double> xTemp, yTemp;
	set<double>::const_iterator it;
	ifstream inpSreda("sreda.txt");
	inpSreda >> numberFields;
	double areaOfSreda = 0;
	for (i = 0; i < numberFields; i++)
	{
		inpSreda >> tempField.x1 >> tempField.x2 >> tempField.y1 >> tempField.y2 >> tempField.lambda >> tempField.gamma;
		sreda.push_back(tempField);
		xTemp.insert(tempField.x1);
		xTemp.insert(tempField.x2);
		yTemp.insert(tempField.y1);
		yTemp.insert(tempField.y2);
		areaOfSreda += (tempField.x2 - tempField.x1) * (tempField.y2 - tempField.y1);
	}
	GenerateNetLikeTelma(xTemp, inpSreda);
	GenerateNetLikeTelma(yTemp, inpSreda);
	set<double>::reverse_iterator rit;
	it = xTemp.begin();
	leftX = *it;
	rit = xTemp.rbegin();
	rightX = *rit;

	it = yTemp.begin();
	leftY = *it;
	rit = yTemp.rbegin();
	rightY = *rit;


	if (fabs(areaOfSreda - (rightX - leftX) * (rightY - leftY)) > 1e-15)
	{
		cout << "Ошибка. Не правильно заданы подобласти. Общая площадь не равна сумме площадей подобластей." << endl;
		system("pause");
		exit(1);
	}

	nX = xTemp.size();
	xNet = new double[nX];
	nY = yTemp.size();
	yNet = new double[nY];
	i = 0;
	for (set<double>::const_iterator it = xTemp.begin(); it != xTemp.end(); it++ , i++)
	{
		xNet[i] = *it;
	}
	i = 0;
	for (set<double>::const_iterator it = yTemp.begin(); it != yTemp.end(); it++ , i++)
	{
		yNet[i] = *it;
	}

	return 0;
}

int FindAreaNumber(int nodes[])
{
	int i;
	for (i = 0; i < sreda.size(); i++)
	{
		if (xy[nodes[0]].x >= sreda[i].x1 && xy[nodes[1]].x <= sreda[i].x2 && xy[nodes[0]].y >= sreda[i].y1 && xy[nodes[2]].y <= sreda[i].y2)
			return i;
	}
	cout << "Ошибка в FindAreaNumber: не найдена подобласть." << endl;
	system("pause");
	exit(1);
	return -1;
}

void outputXYandKE()
{
	int i, j;
	ofstream fileXY("xy.txt");
	ofstream fileNvtr("nvtr.txt");

	//формируем массивы регулярных и терминальных вершин
	for (j = 0; j < nY; j++)
	{
		for (i = 0; i < nX; i++)
		{
			xy.push_back(point(xNet[i], yNet[j]));
		}
	}

	countOfNodes = xy.size();

	//формируем файл xy.txt
	fileXY << xy.size() << endl;
	for (i = 0; i < xy.size(); i++)
	{
		fileXY << xy[i].x << " " << xy[i].y << endl;
	}

	//формируем структуру КЭ
	nvtr tempNvtr;
	for (j = 0; j < nY - 1; j++)
	{
		for (i = 0; i < nX - 1; i++)
		{
			tempNvtr.uzel[0] = NumberNode(xNet[i], yNet[j]);
			tempNvtr.uzel[1] = NumberNode(xNet[i + 1], yNet[j]);
			tempNvtr.uzel[2] = NumberNode(xNet[i], yNet[j + 1]);
			tempNvtr.uzel[3] = NumberNode(xNet[i + 1], yNet[j + 1]);
			tempNvtr.numberField = FindAreaNumber(tempNvtr.uzel);
			KE.push_back(tempNvtr);
		}
	}

	//формируем файл nvtr.txt
	fileNvtr << KE.size() << endl;
	for (i = 0; i < KE.size(); i++)
	{
		for (j = 0; j < 4; j++)
		{
			fileNvtr << KE[i].uzel[j] << " ";
		}
		fileNvtr << endl;
	}
}

void GeneratePortrait()
{
	set<size_t>* portrait = new set<size_t>[countOfNodes];
	for (size_t k = 0; k < KE.size(); k++)
	{
		for (size_t i = 0; i < 4; i++)
		{
			size_t a = KE[k].uzel[i];
			for (size_t j = 0; j < i; j++)
			{
				size_t b = KE[k].uzel[j];
				// Если оба узла не являются терминальными
				if (a < countOfNodes && b < countOfNodes)
				{
					if (b > a)
						portrait[b].insert(a);
					else
						portrait[a].insert(b);
				}
			}
		}
	}

	size_t gg_size = 0;
	for (size_t i = 0; i < countOfNodes; i++)
		gg_size += portrait[i].size();

	countOfGlobalBasicFunc = 4 * countOfNodes;

	ig = new int[countOfGlobalBasicFunc + 1];
	di = new double[countOfGlobalBasicFunc];
	b = new double[countOfGlobalBasicFunc];
	q = new double[countOfGlobalBasicFunc];
	ig[0] = ig[1] = 0;
	for (size_t i = 0; i < countOfGlobalBasicFunc; i++)
	{
		di[i] = b[i] = q[i] = 0;
	}
	for (int i = 0; i < countOfNodes; i++)
	{
		size_t size = portrait[i].size() * 4;
		ig[4 * i + 1] = ig[4 * i] + size;
		ig[4 * i + 2] = ig[4 * i + 1] + size + 1;
		ig[4 * i + 3] = ig[4 * i + 2] + size + 2;
		ig[4 * i + 4] = ig[4 * i + 3] + size + 3;
	}

	jg = new int[ig[countOfGlobalBasicFunc]];
	ggl = new double[ig[countOfGlobalBasicFunc]];
	ggu = new double[ig[countOfGlobalBasicFunc]];

	for (size_t i = 0; i < ig[countOfGlobalBasicFunc]; i++)
	{
		ggl[i] = 0;
		ggu[i] = 0;
	}

	size_t tmp = 0;
	for (size_t i = 0; i < countOfNodes; i++)
	{
		for (int iErm = 4 * i; iErm < 4 * i + 4; iErm++)
		{
			for (set<size_t>::iterator j = portrait[i].begin(); j != portrait[i].end(); ++j)
			{
				for (int jErm = 0; jErm < 4; jErm++)
				{
					jg[tmp] = 4 * *j + jErm;
					tmp++;
				}
			}
			//portrait[i].clear();

			for (int iDop = 4 * i; iDop < iErm; iDop++)
			{
				jg[tmp] = iDop;
				tmp++;
			}
		}
	}
	delete[] portrait;

	int i, j;
	ofstream igOut("ig.txt");
	ofstream jgOut("jg.txt");
	for (i = 0; i <= countOfGlobalBasicFunc; i++)
	{
		igOut << ig[i] << "   " << i << "   ";
		if (i != countOfGlobalBasicFunc)
			igOut << ig[i + 1] - ig[i] << endl;
	}
	cout << endl;
	for (j = 0; j < countOfGlobalBasicFunc; j++)
	{
		jgOut << "string " << j << ": ";
		for (i = ig[j]; i < ig[j + 1]; i++)
		{
			jgOut << jg[i] << " ";
		}
		jgOut << endl;
	}
}

double AnaliticSolution(double x, double y)
{
	//return x;
	//return 5;
	//return x*x + y*y;
	return 1+x+y+x*x+y*y+x*y+x*x*x + y*y*y;
	//return 4 * x * x * x * y * y * y + 5 * x * x * y + 2 * y * y + 10;
}

double AnaliticSolution(point p)
{
	return AnaliticSolution(p.x, p.y);
}

double AnaliticSolutionDx(double x, double y)
{
	double t = (AnaliticSolution(x + hForProizv, y) - AnaliticSolution(x - hForProizv, y)) / (2.0 * hForProizv);
	//double x1 = 12 * x * x * y * y * y + 10 * x * y;
	return t;
}

double AnaliticSolutionDx(point p)
{
	return AnaliticSolutionDx(p.x, p.y);
}

double AnaliticSolutionDy(double x, double y)
{
	double t = (AnaliticSolution(x, y + hForProizv) - AnaliticSolution(x, y - hForProizv)) / (2.0 * hForProizv);
	//double y1 = 12 * x*x*x*y*y + 5 * x*x + 4 * y;
	return t;
}

double AnaliticSolutionDy(point p)
{
	return AnaliticSolutionDy(p.x, p.y);
}

double AnaliticSolutionDxDy(double x, double y)
{
	//return 0;
	return 1;
	//double value = 36 * x * x * y * y + 10 * x;

	//chislennoe
	//double result = (AnaliticSolutionDy(x+hForProizv, y) - AnaliticSolutionDy(x-hForProizv, y)) / (2 * hForProizv);
}

double AnaliticSolutionDxDy(point p)
{
	return AnaliticSolutionDxDy(p.x, p.y);
}

double Lambda(int ielem)
{
	return sreda[KE[ielem].numberField].lambda;
	//return 1;
	//return 5;
}

double Gamma(int ielem)
{
	return sreda[KE[ielem].numberField].gamma;
	//return 1;
	//return 1;
}

double Func(double x,double y)
{
	//return 5;
	//return -4 + x*x + y*y;
	//return -6 * x - 6 * y + x*x*x + y*y*y;
	return -4 -6 * x - 6 * y + 1 + x + y + x*x + y*y+x*y + x*x*x + y*y*y;
	//return -4;
	//return x;
	//return -24 * (x * y * y * y + x * x * x * y) - 10 * y - 4 + 2 * (4 * x * x * x * y * y * y + 5 * x * x * y + 2 * y * y + 10);
}

double Func(int ielem, int localNumber)
{
	double x = xy[KE[ielem].uzel[localNumber]].x;
	double y = xy[KE[ielem].uzel[localNumber]].y;
	return Func(x, y);
}

double FuncDx(int ielem, int localNumber)
{
	double x = xy[KE[ielem].uzel[localNumber]].x;
	double y = xy[KE[ielem].uzel[localNumber]].y;
	double test = -24 * (y * y * y + 3 * x * x * y) + 2 * (12 * x * x * y * y * y + 10 * x * y);
	double chisl = (Func(x + hForProizv, y) - Func(x - hForProizv, y)) / (2.0 * hForProizv);
	return chisl;
}

double FuncDy(int ielem, int localNumber)
{
	double x = xy[KE[ielem].uzel[localNumber]].x;
	double y = xy[KE[ielem].uzel[localNumber]].y;
	double test = -24 * (x * 3 * y * y + x * x * x) - 10 + 2 * (12 * x * x * x * y * y + 5 * x * x + 4 * y);
	double chisl = (Func(x, y + hForProizv) - Func(x, y - hForProizv)) / (2.0 * hForProizv);
	return chisl;
}

double FuncDy(double x,double y)
{
	return (Func(x, y + hForProizv) - Func(x, y - hForProizv)) / (2.0 * hForProizv);
}

double FuncDxDy(int ielem, int localNumber)
{
	double x = xy[KE[ielem].uzel[localNumber]].x;
	double y = xy[KE[ielem].uzel[localNumber]].y;
	double test = -24 * (3 * y * y + 3 * x * x) + 2 * (36 * x * x * y * y + 10 * x);
	double chisl = (FuncDy(x + hForProizv, y) - FuncDy(x - hForProizv, y)) / (2.0 * hForProizv);
	return chisl;
}

double BasicFunc1d(int num, double x, double h)
{
	switch (num)
	{
	case 0: return 1 - 3 * x * x + 2 * x * x * x;
	case 1: return h * (x - 2 * x * x + x * x * x);
	case 2: return 3 * x * x - 2 * x * x * x;
	case 3: return h * (-x * x + x * x * x);
	default:
		{
			cerr << "Error in Basic Func" << endl;
			system("pause");
			exit(1);
		}
	}
}

double** createLocalG1d(double h)
{
	double** loc = new double *[4];
	for (int i = 0; i < 4; i++)
		loc[i] = new double[4];

	double localG[4][4] = {
		{36,3 * h,-36,3 * h},
		{3 * h,4 * h * h,-3 * h,-h * h},
		{-36,-3 * h,36,-3 * h},
		{3 * h,- h * h,-3 * h,4 * h * h}
	};
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
			loc[i][j] = localG[i][j] / (30. * h);
	}
	return loc;
}

double** createLocalM1d(double h)
{
	double** loc = new double *[4];
	for (int i = 0; i < 4; i++)
		loc[i] = new double[4];

	double localM[4][4] = {
		{156,22 * h,54,-13 * h},
		{22 * h,4 * h * h,13 * h,-3 * h * h},
		{54,13 * h,156,-22 * h},
		{-13 * h,-3 * h * h,-22 * h,4 * h * h}
	};
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
			loc[i][j] = localM[i][j] * h / 420.;
	}
	return loc;
}

void CreateLocalMatrixs(int ielem)
{
	double hXlocal, hYlocal;
	locateOfPoint pnt;
	int i, j;
	hXlocal = xy[KE[ielem].uzel[1]].x - xy[KE[ielem].uzel[0]].x;
	hYlocal = xy[KE[ielem].uzel[2]].y - xy[KE[ielem].uzel[0]].y;
	double** localG1dX = createLocalG1d(hXlocal);
	double** localM1dX = createLocalM1d(hXlocal);
	double** localG1dY = createLocalG1d(hYlocal);
	double** localM1dY = createLocalM1d(hYlocal);

	for (i = 0; i < 16; i++)
	{
		localB[i] = 0;
		int muI = 2 * ((i / 4) % 2) + (i % 2);
		int nuI = 2 * (i / 8) + ((i / 2) % 2);
		for (j = 0; j < 16; j++)
		{
			int muJ = 2 * ((j / 4) % 2) + (j % 2);
			int nuJ = 2 * (j / 8) + ((j / 2) % 2);

			localMatrix[i][j] = Lambda(ielem) * (localG1dX[muI][muJ] * localM1dY[nuI][nuJ] + localM1dX[muI][muJ] * localG1dY[nuI][nuJ]);
			localMatrix[i][j] += Gamma(ielem) * localM1dX[muI][muJ] * localM1dY[nuI][nuJ];

			switch (j % 4)
			{
			case 0: localB[i] += localM1dX[muI][muJ] * localM1dY[nuI][nuJ] * Func(ielem, j / 4);
				break;
			case 1: localB[i] += localM1dX[muI][muJ] * localM1dY[nuI][nuJ] * FuncDx(ielem, j / 4);
				break;
			case 2: localB[i] += localM1dX[muI][muJ] * localM1dY[nuI][nuJ] * FuncDy(ielem, j / 4);
				break;
			case 3: localB[i] += localM1dX[muI][muJ] * localM1dY[nuI][nuJ] * FuncDxDy(ielem, j / 4);
				break;
			}
		}
	}
}

void AddToMatrix(int posI, int posJ, double el)
{
	int tmp;
	if (posI == posJ)
	{
		di[posI] += el;
		return;
	}
	else
	{
		if (posI < posJ)
		{
			return;
			tmp = posI;
			posI = posJ;
			posJ = tmp;
		}
		for (tmp = ig[posI]; tmp < ig[posI + 1]; tmp++)
		{
			if (jg[tmp] == posJ)
			{
				ggl[tmp] += el;
				return;
			}
		}
	}
}

void Addition(int ielem)
{
	int i, j, iSchet, jSchet;
	for (i = 0; i < 4; i++)
	{
		for (iSchet = 0; iSchet < 4; iSchet++)
		{
			b[4 * KE[ielem].uzel[i] + iSchet] += localB[4 * i + iSchet];
			for (j = 0; j <= i; j++)
			{
				for (jSchet = 0; jSchet < 4; jSchet++)
				{
					AddToMatrix(4 * KE[ielem].uzel[i] + iSchet, 4 * KE[ielem].uzel[j] + jSchet, localMatrix[4 * i + iSchet][4 * j + jSchet]);
				}
			}
		}
	}
}

void PrintPlotMatrix(bool flag_simmeric)
{
	int i, j;
	double** APlot = new double*[countOfGlobalBasicFunc];
	for (i = 0; i < countOfGlobalBasicFunc; i++)
	{
		APlot[i] = new double[countOfGlobalBasicFunc];
		for (j = 0; j < countOfGlobalBasicFunc; j++)
		{
			APlot[i][j] = 0;
		}
	}
	if (flag_simmeric)
		for (i = 0; i < countOfGlobalBasicFunc; i++)
		{
			APlot[i][i] = di[i];
			for (j = ig[i]; j < ig[i + 1]; j++)
			{
				APlot[i][jg[j]] = ggl[j];
				APlot[jg[j]][i] = ggl[j];
			}
		}
	else
		for (i = 0; i < countOfGlobalBasicFunc; i++)
		{
			APlot[i][i] = di[i];
			for (j = ig[i]; j < ig[i + 1]; j++)
			{
				APlot[i][jg[j]] = ggl[j];
				APlot[jg[j]][i] = ggu[j];
				//APlot[jg[j]][i] = 0;
			}
		}

	for (i = 0; i < countOfGlobalBasicFunc; i++)
	{
		for (j = 0; j < countOfGlobalBasicFunc; j++)
		{
			outFile << setw(15) << APlot[i][j];
		}
		outFile << endl;
	}
	outFile << endl;

	for (i = 0; i < countOfGlobalBasicFunc; i++)
	{
		outFile << setw(15) << b[i];
	}
	outFile << endl;
	outFile << endl;
}


//intXorY can be 0 or 1
void doEdge1(ofstream& outEdge1File, int intXorY, int kolvoRegularNode, double* varNet1, int nVarNet1, double uncnownValuePlosk)
{
	for (int iVar1 = 0; iVar1 < nVarNet1; iVar1++)
	{
		point goal(0, 0);
		int iProizv = 0;
		double valueProizv = 0;
		if (intXorY == 0)
		{
			goal = point(uncnownValuePlosk, varNet1[iVar1]);
			iProizv = 2;
			valueProizv = AnaliticSolutionDy(goal);
		}
		else
		{
			goal = point(varNet1[iVar1], uncnownValuePlosk);
			iProizv = 1;
			valueProizv = AnaliticSolutionDx(goal);
		}
		int k = 4 * indexXY(goal);
		di[k] = di[k + iProizv] = 1;
		b[k] = AnaliticSolution(goal);
		b[k + iProizv] = valueProizv;
		for (int iter = 0, prirash = 0; iter < 2; iter++)
		{
			for (int m = ig[k + prirash]; m < ig[k + prirash + 1]; m++)
			{
				ggl[m] = 0;
			}
			prirash = iProizv;
		}
		for (int l = 0; l < kolvoRegularNode; l++)
		{
			for (int m = ig[l]; m < ig[l + 1]; m++)
			{
				if (k == jg[m] || jg[m] == k + iProizv)
				{
					ggu[m] = 0;
				}
			}
		}
		outEdge1File << k << '\t' << b[k] << endl;
		outEdge1File << k + iProizv << '\t' << b[k + iProizv] << endl;
	}
}

//
void Edge1_not_sim(bool up, bool down, bool left, bool right)
{
	ofstream ku1("ku1.txt");
	if (down)
		doEdge1(ku1, 1, countOfGlobalBasicFunc, xNet, nX, yNet[0]);
	if (up)
		doEdge1(ku1, 1, countOfGlobalBasicFunc, xNet, nX, yNet[nY - 1]);
	if (left)
		doEdge1(ku1, 0, countOfGlobalBasicFunc, yNet, nY, xNet[0]);
	if (right)
		doEdge1(ku1, 0, countOfGlobalBasicFunc, yNet, nY, xNet[nX - 1]);
}

//intXorY can be 0 or 1
//0 - если не меняется Х, иначе 1
void doEdge2(ofstream& outEdge2File, int intXorY, int normalDirect, double* varNet1, int nVarNet1, double uncnownValuePlosk)
{
	double dUdn[4], dh1;
	int indNodes[4];
	double** M1d;
	for (int iVar1 = 0; iVar1 < nVarNet1 - 1; iVar1++)
	{
		point goal(0, 0);
		if (intXorY == 0)
		{
			goal = point(uncnownValuePlosk, varNet1[iVar1]);
			point nextNode(uncnownValuePlosk, varNet1[iVar1 + 1]);
			dUdn[0] = normalDirect * AnaliticSolutionDx(goal);
			dUdn[1] = normalDirect * AnaliticSolutionDxDy(goal);
			dUdn[2] = normalDirect * AnaliticSolutionDx(nextNode);
			dUdn[3] = normalDirect * AnaliticSolutionDxDy(nextNode);
			indNodes[0] = 4 * indexXY(goal);
			indNodes[1] = 4 * indexXY(goal) + 2;
			indNodes[2] = 4 * indexXY(nextNode);
			indNodes[3] = 4 * indexXY(nextNode) + 2;
		}
		else
		{
			goal = point(varNet1[iVar1], uncnownValuePlosk);
			point nextNode(varNet1[iVar1 + 1], uncnownValuePlosk);
			dUdn[0] = normalDirect * AnaliticSolutionDy(goal);
			dUdn[1] = normalDirect * AnaliticSolutionDxDy(goal);
			dUdn[2] = normalDirect * AnaliticSolutionDy(nextNode);
			dUdn[3] = normalDirect * AnaliticSolutionDxDy(nextNode);
			indNodes[0] = 4 * indexXY(goal);
			indNodes[1] = 4 * indexXY(goal) + 1;
			indNodes[2] = 4 * indexXY(nextNode);
			indNodes[3] = 4 * indexXY(nextNode) + 1;
		}
		int indKE = FindNumberOfKe(indNodes[0]/4, indNodes[2]/4);
		dh1 = varNet1[iVar1 + 1] - varNet1[iVar1];
		M1d = createLocalM1d(dh1);

		for (size_t iter1 = 0; iter1 < 4; iter1++)
		{
			double value = 0;
			for (size_t iter2 = 0; iter2 < 4; iter2++)
			{
				value += dUdn[iter2] * M1d[iter1][iter2] * Lambda(indKE);
			}
			b[indNodes[iter1]] += value;
			delete[]M1d[iter1];
			outEdge2File << setw(10) << indNodes[iter1] << setw(15) << value << endl;
		}
		delete[]M1d;
	}
}

void Edge2_not_sim(bool up, bool down, bool left, bool right)
{
	ofstream ku2("ku2.txt");
	if (down)
		doEdge2(ku2, 1, -1, xNet, nX, yNet[0]);
	if (up)
		doEdge2(ku2, 1, 1, xNet, nX, yNet[nY - 1]);
	if (left)
		doEdge2(ku2, 0, -1, yNet, nY, xNet[0]);
	if (right)
		doEdge2(ku2, 0, 1, yNet, nY, xNet[nX - 1]);
}

//intXorY can be 0 or 1
//0 - если не меняется Х, иначе 1
void doEdge3(ofstream& outEdge3File, int intXorY, int normalDirect, double* varNet1, int nVarNet1, double uncnownValuePlosk)
{
	double dUdn[4], dh1, ubetta[4];
	int indNodes[4];
	double** M1d;
	for (int iVar1 = 0; iVar1 < nVarNet1 - 1; iVar1++)
	{
		point goal(0, 0);
		int indKE;
		if (intXorY == 0)
		{
			goal = point(uncnownValuePlosk, varNet1[iVar1]);
			point nextNode(uncnownValuePlosk, varNet1[iVar1 + 1]);
			dUdn[0] = normalDirect * AnaliticSolutionDx(goal);
			dUdn[1] = normalDirect * AnaliticSolutionDxDy(goal);
			dUdn[2] = normalDirect * AnaliticSolutionDx(nextNode);
			dUdn[3] = normalDirect * AnaliticSolutionDxDy(nextNode);
			indNodes[0] = 4 * indexXY(goal);
			indNodes[1] = 4 * indexXY(goal) + 2;
			indNodes[2] = 4 * indexXY(nextNode);
			indNodes[3] = 4 * indexXY(nextNode) + 2;

			indKE = FindNumberOfKe(indNodes[0] / 4, indNodes[2] / 4);

			double a1 = normalDirect * AnaliticSolutionDxDy(goal.x, goal.y + hForProizv)* Lambda(indKE) + AnaliticSolution(goal.x, goal.y + hForProizv);
			double a2 = (normalDirect * AnaliticSolutionDxDy(goal.x, goal.y - hForProizv)* Lambda(indKE) + AnaliticSolution(goal.x, goal.y - hForProizv));
			ubetta[1] = (a1 - a2) / (2 * hForProizv);
			ubetta[3] = (normalDirect * AnaliticSolutionDxDy(nextNode.x, nextNode.y + hForProizv)* Lambda(indKE) + AnaliticSolution(nextNode.x, nextNode.y + hForProizv) -
				(normalDirect * AnaliticSolutionDxDy(nextNode.x, nextNode.y - hForProizv)* Lambda(indKE) + AnaliticSolution(nextNode.x, nextNode.y - hForProizv))) / (2 * hForProizv);
		}
		else
		{
			goal = point(varNet1[iVar1], uncnownValuePlosk);
			point nextNode(varNet1[iVar1 + 1], uncnownValuePlosk);
			dUdn[0] = normalDirect * AnaliticSolutionDy(goal);
			dUdn[1] = normalDirect * AnaliticSolutionDxDy(goal);
			dUdn[2] = normalDirect * AnaliticSolutionDy(nextNode);
			dUdn[3] = normalDirect * AnaliticSolutionDxDy(nextNode);
			indNodes[0] = 4 * indexXY(goal);
			indNodes[1] = 4 * indexXY(goal) + 1;
			indNodes[2] = 4 * indexXY(nextNode);
			indNodes[3] = 4 * indexXY(nextNode) + 1;

			indKE = FindNumberOfKe(indNodes[0] / 4, indNodes[2] / 4);

			ubetta[1] = (normalDirect * AnaliticSolutionDxDy(goal.x + hForProizv, goal.y)* Lambda(indKE) + AnaliticSolution(goal.x + hForProizv, goal.y) -
				(normalDirect * AnaliticSolutionDxDy(goal.x - hForProizv, goal.y)* Lambda(indKE) + AnaliticSolution(goal.x - hForProizv, goal.y))) / (2 * hForProizv);
			ubetta[3] = (normalDirect * AnaliticSolutionDxDy(nextNode.x + hForProizv, nextNode.y)* Lambda(indKE) + AnaliticSolution(nextNode.x + hForProizv, nextNode.y) -
				(normalDirect * AnaliticSolutionDxDy(nextNode.x - hForProizv, nextNode.y)* Lambda(indKE) + AnaliticSolution(nextNode.x - hForProizv, nextNode.y))) / (2 * hForProizv);
		}
		dh1 = varNet1[iVar1 + 1] - varNet1[iVar1];
		M1d = createLocalM1d(dh1);
		double betta = 1;

		for (size_t i = 0; i < 4; i+=2)
		{
			ubetta[i] = AnaliticSolution(xy[indNodes[i]/4]) + dUdn[i];
		}

		for (size_t iter1 = 0; iter1 < 4; iter1++)
		{
			double value = 0;
			for (size_t iter2 = 0; iter2 < 4; iter2++)
			{
				value += betta*ubetta[iter2] * M1d[iter1][iter2];
				AddToMatrix(indNodes[iter1], indNodes[iter2], betta*M1d[iter1][iter2]);
			}
			b[indNodes[iter1]] += value;
			delete[]M1d[iter1];
			outEdge3File << setw(10) << indNodes[iter1] << setw(15) << value << endl;
		}
		delete[]M1d;
	}
}

void Edge3_not_sim(bool up, bool down, bool left, bool right)
{
	ofstream ku3("ku3.txt");
	if (down)
		doEdge3(ku3, 1, -1, xNet, nX, yNet[0]);
	if (up)
		doEdge3(ku3, 1, 1, xNet, nX, yNet[nY - 1]);
	if (left)
		doEdge3(ku3, 0, -1, yNet, nY, xNet[0]);
	if (right)
		doEdge3(ku3, 0, 1, yNet, nY, xNet[nX - 1]);
}

void PrintLocalMatrix()
{
	int i, j;
	for (i = 0; i < 4; i++)
	{
		for (j = 0; j < 4; j++)
		{
			outFile << setw(15) << localMatrix[i][j];
		}
		outFile << endl;
	}
	outFile << endl;
}

void TestFromBook()
{
	double** Mx = createLocalM1d(10);
	double** My = createLocalM1d(xy[3].y - xy[0].y);

	{
		double vect[4] = {-500,-100,-2000,-200};
		double bVect[4];
		int konvertation[4] = {0,1,4,5};

		//второе краевое S12
		for (int i = 0; i < 4; i++)
		{
			bVect[i] = 0;
			for (int j = 0; j < 4; j++)
			{
				bVect[i] += Mx[i][j] * vect[j];
			}
			b[konvertation[i]] += bVect[i];
		}
	}

	//второе краевое S22
	{
		double vect[4] = {0,-100,-11.2,-136};
		double bVect[4];
		int konvertation[4] = {0,2,8,10};
		for (int i = 0; i < 4; i++)
		{
			bVect[i] = 0;
			for (int j = 0; j < 4; j++)
			{
				bVect[i] += My[i][j] * vect[j];
			}
			b[konvertation[i]] += bVect[i];
		}
	}

	//третье краевое S3
	{
		double vect[4] = {10,2200,266.82,3304.4};
		double bVect[4], A[4][4];
		int konvertation[4] = {4,6,12,14};
		for (int i = 0; i < 4; i++)
		{
			bVect[i] = 0;
			for (int j = 0; j < 4; j++)
			{
				bVect[i] += My[i][j] * vect[j];
				A[i][j] = My[i][j];
				AddToMatrix(konvertation[i], konvertation[j], A[i][j]);
				//AddToMatrix(konvertation[j], konvertation[i], A[i][j]);
			}
			b[konvertation[i]] += bVect[i];
		}
	}

	//первое краевое
	{
		int konvertation[4] = {8,9,12,13};
		for (int i = 0; i < 4; i++)
		{
			di[konvertation[i]] = 1;
			if (konvertation[i] == 8 || konvertation[i] == 12)
				b[konvertation[i]] = AnaliticSolution(xNet[i / 2], yNet[1]);
			else
				b[konvertation[i]] = AnaliticSolutionDx(xNet[i / 2], yNet[1]);
			for (int m = ig[konvertation[i]]; m < ig[konvertation[i] + 1]; m++)
			{
				ggl[m] = 0;
			}
			for (int l = 0; l < countOfGlobalBasicFunc; l++)
			{
				for (int m = ig[l]; m < ig[l + 1]; m++)
				{
					if (konvertation[i] == jg[m])
					{
						ggu[m] = 0;
					}
				}
			}
		}
	}
}

void GenerateMatrix()
{
	int ielem, i, j;
	for (ielem = 0; ielem < KE.size(); ielem++)
	{
		CreateLocalMatrixs(ielem);
		Addition(ielem);
		//PrintLocalMatrix();
		//PrintPlotMatrix(true);
	}
	//PrintPlotMatrix(false);

	//Edge2_not_sim(1, 1, 1, 1);
	for (i = 0; i < ig[countOfGlobalBasicFunc]; i++)
	{
		ggu[i] = ggl[i];
	}
	Edge1_not_sim(1, 1, 1, 1);
	//TestFromBook();
	PrintPlotMatrix(false);
	//Edge1_not_sim(1, 1, 1, 1);
	//PrintPlotMatrix(true);
}

void mult(double* res, double* v)
{
	for (int i = 0; i < countOfNodes; i++)
		res[i] = 0;
	for (int i = 0; i < countOfNodes; i++)
	{
		for (int j = ig[i]; j < ig[i + 1]; j++)
		{
			res[i] += ggl[j] * v[jg[j]];
			res[jg[j]] += ggu[j] * v[i];
		}
		res[i] += di[i] * v[i];
	}
}

double ScalarMult(double* v1, double* v2, int dimQ)
{
	int i;
	double result;
	result = 0;
	for (i = 0; i < dimQ; i++)
	{
		result += v1[i] * v2[i];
	}
	return result;
}

//	Умножение матрицы на вектор
void MultMatrixOnVector(double* in, double* out, int dimQ)
{
	int i, j;
	double* out1;
	out1 = new double[dimQ];
	for (i = 0; i < dimQ; i++)
	{
		out1[i] = di[i] * in[i];
		for (j = ig[i]; j < ig[i + 1]; j++)
		{
			out1[i] += ggl[j] * in[jg[j]];
			out1[jg[j]] += ggu[j] * in[i];
		}
	}
	for (i = 0; i < dimQ; i++)
		out[i] = out1[i];
	delete[] out1;
}

void runMSG(int dimQ)
{
	int maxiter = 10000, i;
	double alfa, alfachisl, alfaznam, beta, betachisl, betaznam, checkE, epsMSG = 1e-16;
	double* r = new double[dimQ];
	double* s = new double[dimQ];
	double* z = new double[dimQ];
	double* rout = new double[dimQ];
	for (i = 0; i < dimQ; i++)
	{
		s[i] = rout[i] = r[i] = q[i] = z[i] = 0;
	}
	MultMatrixOnVector(q, r, dimQ);
	for (i = 0; i < dimQ; i++)
	{
		r[i] = b[i] - r[i];
		z[i] = r[i];
	}
	checkE = sqrt(ScalarMult(r, r, dimQ) / ScalarMult(b, b, dimQ));
	for (iter = 0; iter < maxiter && checkE >= epsMSG; iter++)
	{
		alfachisl = ScalarMult(r, r, dimQ);
		MultMatrixOnVector(z, s, dimQ);
		alfaznam = ScalarMult(s, z, dimQ);
		alfa = alfachisl / alfaznam;
		for (i = 0; i < dimQ; i++)
		{
			q[i] = q[i] + alfa * z[i];
			rout[i] = r[i] - alfa * s[i];
		}
		betachisl = ScalarMult(rout, rout, dimQ);
		betaznam = ScalarMult(r, r, dimQ);
		beta = betachisl / betaznam;
		for (i = 0; i < dimQ; i++)
		{
			z[i] = rout[i] + beta * z[i];
			r[i] = rout[i];
		}
		checkE = sqrt(ScalarMult(r, r, dimQ) / ScalarMult(b, b, dimQ));
	}
}

void runLOS(int dimQ)
{
	int maxiter = 10000, i;
	double alfa, alfachisl, alfaznam, beta, betachisl, betaznam, checkE, epsMSG = 1e-16, startNeviazka;
	double* r = new double[dimQ];
	double* s = new double[dimQ];
	double* z = new double[dimQ];
	double* p = new double[dimQ];
	double* rout = new double[dimQ];
	for (i = 0; i < dimQ; i++)
	{
		s[i] = rout[i] = r[i] = q[i] = z[i] = p[i] = 0;
	}
	MultMatrixOnVector(q, r, dimQ);
	for (i = 0; i < dimQ; i++)
	{
		r[i] = b[i] - r[i];
		z[i] = r[i];
	}
	MultMatrixOnVector(z, p, dimQ);
	checkE = sqrt(ScalarMult(r, r, dimQ) / ScalarMult(b, b, dimQ));
	//startNeviazka = checkE = ScalarMult(r, r);
	for (iter = 0; iter < maxiter && checkE >= epsMSG; iter++)
	{
		alfachisl = ScalarMult(p, r, dimQ);
		alfaznam = ScalarMult(p, p, dimQ);
		alfa = alfachisl / alfaznam;
		for (i = 0; i < dimQ; i++)
		{
			q[i] = q[i] + alfa * z[i];
			r[i] = r[i] - alfa * p[i];
		}
		MultMatrixOnVector(r, rout, dimQ);
		betachisl = ScalarMult(p, rout, dimQ);
		betaznam = ScalarMult(p, p, dimQ);
		beta = -betachisl / betaznam;
		for (i = 0; i < dimQ; i++)
		{
			z[i] = r[i] + beta * z[i];
			p[i] = rout[i] + beta * p[i];
		}
		checkE = sqrt(ScalarMult(r, r, dimQ) / ScalarMult(b, b, dimQ));
	}
}

void outputResult()
{
	int i;
	double otnosPogreshnost = 0, chislitPogreshnosti = 0, analitSolvInNode;
	outFile << "iter = " << iter << endl;
	for (i = 0; i < countOfNodes; i++)
	{
		analitSolvInNode = AnaliticSolution(xy[i].x, xy[i].y);
		outFile << setw(15) << q[4 * i] << "\t" << setw(15) << analitSolvInNode << "\t" << setw(15) << abs(q[4 * i] - analitSolvInNode) << endl;
		chislitPogreshnosti += (q[4 * i] - analitSolvInNode) * (q[4 * i] - analitSolvInNode);
		otnosPogreshnost += analitSolvInNode * analitSolvInNode;

		analitSolvInNode = AnaliticSolutionDx(xy[i].x, xy[i].y);
		outFile << setw(15) << q[4 * i + 1] << "\t" << setw(15) << analitSolvInNode << "\t" << setw(15) << abs(q[4 * i + 1] - analitSolvInNode) << endl;
		chislitPogreshnosti += (q[4 * i + 1] - analitSolvInNode) * (q[4 * i + 1] - analitSolvInNode);
		otnosPogreshnost += analitSolvInNode * analitSolvInNode;

		analitSolvInNode = AnaliticSolutionDy(xy[i].x, xy[i].y);
		outFile << setw(15) << q[4 * i + 2] << "\t" << setw(15) << analitSolvInNode << "\t" << setw(15) << abs(q[4 * i + 2] - analitSolvInNode) << endl;
		chislitPogreshnosti += (q[4 * i + 2] - analitSolvInNode) * (q[4 * i + 2] - analitSolvInNode);
		otnosPogreshnost += analitSolvInNode * analitSolvInNode;

		analitSolvInNode = AnaliticSolutionDxDy(xy[i].x, xy[i].y);
		outFile << setw(15) << q[4 * i + 3] << "\t" << setw(15) << analitSolvInNode << "\t" << setw(15) << abs(q[4 * i + 3] - analitSolvInNode) << endl;
		chislitPogreshnosti += (q[4 * i + 3] - analitSolvInNode) * (q[4 * i + 3] - analitSolvInNode);
		otnosPogreshnost += analitSolvInNode * analitSolvInNode;
	}
	otnosPogreshnost = sqrt(chislitPogreshnosti / otnosPogreshnost);
	outFile << "Относительная погрешность: " << otnosPogreshnost << endl << "Количество КЭ: " << KE.size() << endl << "Количество узлов: " << countOfNodes;
}

int main(int argc, char* argv[])
{
	setlocale(LC_ALL, "rus");
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGB);

	GenerateNet();

	outputXYandKE();

	GeneratePortrait();
	GenerateMatrix();
	//runMSG();
	runLOS(countOfGlobalBasicFunc);

	outputResult();

	Width = 600;
	Height = 600;
	glutInitWindowSize(Width, Height);
	glutCreateWindow("GLUT");
	glutDisplayFunc(Display);
	glutReshapeFunc(Reshape);
	glutIdleFunc(Display);
	glutMainLoop();
	return 1;
}
