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

//matrix
int *ig, *jg;
double *ggl, *ggu, *di, *b, *q;
int iter;

double localB[4], localMatrix[4][4];
double helpG1[4][4] = {{2, -2, 1, -1},{-2, 2, -1, 1},{1, -1, 2, -2},{-1, 1, -2, 2}};
double helpG2[4][4] = {{2, 1, -2, -1},{1, 2, -1, -2},{-2, -1, 2, 1},{-1, -2, 1, 2}};
double helpM[4][4] = {{4, 2, 2, 1},{2, 4, 1, 2},{2, 1, 4, 2},{1, 2, 2, 4}};

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
	cout << "Ошибка краевых условий. Не найден КЭ." << endl;
	exit(1);
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

	ig = new int[countOfNodes + 1];
	di = new double[countOfNodes];
	b = new double[countOfNodes];
	q = new double[countOfNodes];
	ig[0] = ig[1] = 0;
	for (size_t i = 0; i < countOfNodes; i++)
	{
		di[i] = b[i] = q[i] = 0;
		ig[i + 1] = ig[i] + portrait[i].size();
	}
	jg = new int[ig[countOfNodes]];
	ggl = new double[ig[countOfNodes]];
	ggu = new double[ig[countOfNodes]];

	for (size_t i = 0; i < ig[countOfNodes]; i++)
	{
		ggl[i] = 0;
		ggu[i] = 0;
	}

	size_t tmp = 0;
	for (size_t i = 0; i < countOfNodes; i++)
	{
		for (set<size_t>::iterator j = portrait[i].begin(); j != portrait[i].end(); ++j)
		{
			jg[tmp] = *j;
			tmp++;
		}
		portrait[i].clear();
	}
	delete[] portrait;

	int i, j;
	ofstream igOut("ig.txt");
	ofstream jgOut("jg.txt");
	for (i = 0; i <= countOfNodes; i++)
	{
		igOut << ig[i] << "   " << i << "   ";
		if (i != countOfNodes)
			igOut << ig[i + 1] - ig[i] << endl;
	}
	cout << endl;
	for (j = 0; j < countOfNodes; j++)
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
	return x + y;
	//return x*x + y*y;
}

double Lambda(int ielem)
{
	//return sreda[KE[ielem].numberField].lambda;
	return 1;
	//return 5;
}

double Gamma(int ielem)
{
	//return sreda[KE[ielem].numberField].gamma;
	//return 1;
	return 1;
}

double Func(int ielem, int localNumber)
{
	return xy[KE[ielem].uzel[localNumber]].x + xy[KE[ielem].uzel[localNumber]].y;
	//return sreda[KE[ielem].numberField].gamma*AnaliticSolution(xy[KE[ielem].uzel[localNumber]].x, xy[KE[ielem].uzel[localNumber]].y);
	//return xy[KE[ielem].uzel[localNumber]].x;
	//return sreda[KE[ielem].numberField].gamma*(xy[KE[ielem].uzel[localNumber]].x + xy[KE[ielem].uzel[localNumber]].y + 5);
	//return 1;
	//return -4;
	//return (xy[KE[ielem].uzel[localNumber]].x + xy[KE[ielem].uzel[localNumber]].y) + 5;
	//return 0;
	//return 2*(xy[KE[ielem].uzel[localNumber]].x*xy[KE[ielem].uzel[localNumber]].x + xy[KE[ielem].uzel[localNumber]].y*xy[KE[ielem].uzel[localNumber]].y) - 20;
}

void CreateLocalMatrixs(int ielem)
{
	double hXlocal, hYlocal;
	locateOfPoint pnt;
	int i, j;
	hXlocal = xy[KE[ielem].uzel[1]].x - xy[KE[ielem].uzel[0]].x;
	hYlocal = xy[KE[ielem].uzel[2]].y - xy[KE[ielem].uzel[0]].y;
	for (i = 0; i < 4; i++)
	{
		localB[i] = 0;
		for (j = 0; j < 4; j++)
		{
			localB[i] += hXlocal * hYlocal * helpM[i][j] * Func(ielem, j) / 36.0;
			localMatrix[i][j] = Lambda(ielem) / 6.0 * (hYlocal * helpG1[i][j] / hXlocal + hXlocal * helpG2[i][j] / hYlocal) + Gamma(ielem) * hXlocal * hYlocal / 36.0 * helpM[i][j];
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
	int i, j;
	for (i = 0; i < 4; i++)
	{
		b[KE[ielem].uzel[i]] += localB[i];

		for (j = 0; j < 4; j++)
		{
			AddToMatrix(KE[ielem].uzel[i], KE[ielem].uzel[j], localMatrix[i][j]);
		}
	}
}

void PrintPlotMatrix(bool flag_simmeric)
{
	int i, j;
	double** APlot = new double*[countOfNodes];
	for (i = 0; i < countOfNodes; i++)
	{
		APlot[i] = new double[countOfNodes];
		for (j = 0; j < countOfNodes; j++)
		{
			APlot[i][j] = 0;
		}
	}
	if (flag_simmeric)
		for (i = 0; i < countOfNodes; i++)
		{
			APlot[i][i] = di[i];
			for (j = ig[i]; j < ig[i + 1]; j++)
			{
				APlot[i][jg[j]] = ggl[j];
				APlot[jg[j]][i] = ggl[j];
			}
		}
	else
		for (i = 0; i < countOfNodes; i++)
		{
			APlot[i][i] = di[i];
			for (j = ig[i]; j < ig[i + 1]; j++)
			{
				APlot[i][jg[j]] = ggl[j];
				APlot[jg[j]][i] = ggu[j];
			}
		}

	for (i = 0; i < countOfNodes; i++)
	{
		for (j = 0; j < countOfNodes; j++)
		{
			outFile << setw(15) << APlot[i][j];
		}
		outFile << endl;
	}
	outFile << endl;

	for (i = 0; i < countOfNodes; i++)
	{
		outFile << setw(15) << b[i];
	}
	outFile << endl;
	outFile << endl;
}

void Edge1_sim(bool up, bool down, bool left, bool right)
{
	ofstream ku1("ku1.txt");
	int i, j, k, m, l;

	//--------
	for (i = 0; i < nX; i++)
	{
		if (down == true)
		{
			j = 0;
			k = NumberNode(xNet[i], yNet[j]);
			for (m = ig[k]; m < ig[k + 1]; m++)
			{
				b[jg[m]] -= ggl[m] * AnaliticSolution(xNet[i], yNet[j]);
				ggl[m] = 0;
			}
		}
		if (up == true)
		{
			j = nY - 1;
			k = NumberNode(xNet[i], yNet[j]);
			for (m = ig[k]; m < ig[k + 1]; m++)
			{
				b[jg[m]] -= ggl[m] * AnaliticSolution(xNet[i], yNet[j]);
				ggl[m] = 0;
			}
		}
	}

	for (j = 0; j < nY; j++)
	{
		if (left == true)
		{
			i = 0;
			k = NumberNode(xNet[i], yNet[j]);
			for (m = ig[k]; m < ig[k + 1]; m++)
			{
				b[jg[m]] -= ggl[m] * AnaliticSolution(xNet[i], yNet[j]);
				ggl[m] = 0;
			}
		}
		if (right == true)
		{
			i = nX - 1;
			k = NumberNode(xNet[i], yNet[j]);
			for (m = ig[k]; m < ig[k + 1]; m++)
			{
				b[jg[m]] -= ggl[m] * AnaliticSolution(xNet[i], yNet[j]);
				ggl[m] = 0;
			}
		}
	}

	//------------
	for (i = 0; i < nX; i++)
	{
		if (down == true)
		{
			j = 0;
			k = NumberNode(xNet[i], yNet[j]);
			di[k] = 1;
			b[k] = AnaliticSolution(xNet[i], yNet[j]);
			for (m = ig[k]; m < ig[k + 1]; m++)
			{
				ggl[m] = 0;
			}
			for (l = 0; l < countOfNodes; l++)
			{
				for (m = ig[l]; m < ig[l + 1]; m++)
				{
					if (k == jg[m])
					{
						b[l] -= b[k] * ggl[m];
						ggl[m] = 0;
					}
				}
			}
			ku1 << k << '\t' << b[k] << endl;
		}
		if (up == true)
		{
			j = nY - 1;
			k = NumberNode(xNet[i], yNet[j]);
			di[k] = 1;
			b[k] = AnaliticSolution(xNet[i], yNet[j]);
			for (m = ig[k]; m < ig[k + 1]; m++)
			{
				ggl[m] = 0;
			}
			for (l = 0; l < countOfNodes; l++)
			{
				for (m = ig[l]; m < ig[l + 1]; m++)
				{
					if (k == jg[m])
					{
						b[l] -= b[k] * ggl[m];
						ggl[m] = 0;
					}
				}
			}
			ku1 << k << '\t' << b[k] << endl;
		}
	}

	for (j = 0; j < nY; j++)
	{
		if (left == true)
		{
			i = 0;
			k = NumberNode(xNet[i], yNet[j]);
			di[k] = 1;
			b[k] = AnaliticSolution(xNet[i], yNet[j]);
			for (m = ig[k]; m < ig[k + 1]; m++)
			{
				ggl[m] = 0;
			}
			for (l = 0; l < countOfNodes; l++)
			{
				for (m = ig[l]; m < ig[l + 1]; m++)
				{
					if (k == jg[m])
					{
						b[l] -= b[k] * ggl[m];
						ggl[m] = 0;
					}
				}
			}
			ku1 << k << '\t' << b[k] << endl;
		}
		if (right == true)
		{
			i = nX - 1;
			k = NumberNode(xNet[i], yNet[j]);
			di[k] = 1;
			b[k] = AnaliticSolution(xNet[i], yNet[j]);
			for (m = ig[k]; m < ig[k + 1]; m++)
			{
				ggl[m] = 0;
			}
			for (l = 0; l < countOfNodes; l++)
			{
				for (m = ig[l]; m < ig[l + 1]; m++)
				{
					if (k == jg[m])
					{
						b[l] -= b[k] * ggl[m];
						ggl[m] = 0;
					}
				}
			}
			ku1 << k << '\t' << b[k] << endl;
		}
	}
}

void Edge1_not_sim(bool up, bool down, bool left, bool right)
{
	ofstream ku1("ku1.txt");
	int i, j, k, m, l;

	for (i = 0; i < nX; i++)
	{
		if (down == true)
		{
			j = 0;
			k = NumberNode(xNet[i], yNet[j]);
			di[k] = 1;
			b[k] = AnaliticSolution(xNet[i], yNet[j]);
			for (m = ig[k]; m < ig[k + 1]; m++)
			{
				ggl[m] = 0;
			}
			for (l = 0; l < countOfNodes; l++)
			{
				for (m = ig[l]; m < ig[l + 1]; m++)
				{
					if (k == jg[m])
					{
						ggu[m] = 0;
					}
				}
			}
			ku1 << k << '\t' << b[k] << endl;
		}
		if (up == true)
		{
			j = nY - 1;
			k = NumberNode(xNet[i], yNet[j]);
			di[k] = 1;
			b[k] = AnaliticSolution(xNet[i], yNet[j]);
			for (m = ig[k]; m < ig[k + 1]; m++)
			{
				ggl[m] = 0;
			}
			for (l = 0; l < countOfNodes; l++)
			{
				for (m = ig[l]; m < ig[l + 1]; m++)
				{
					if (k == jg[m])
					{
						ggu[m] = 0;
					}
				}
			}
			ku1 << k << '\t' << b[k] << endl;
		}
	}

	for (j = 0; j < nY; j++)
	{
		if (left == true)
		{
			i = 0;
			k = NumberNode(xNet[i], yNet[j]);
			di[k] = 1;
			b[k] = AnaliticSolution(xNet[i], yNet[j]);
			for (m = ig[k]; m < ig[k + 1]; m++)
			{
				ggl[m] = 0;
			}
			for (l = 0; l < countOfNodes; l++)
			{
				for (m = ig[l]; m < ig[l + 1]; m++)
				{
					if (k == jg[m])
					{
						ggu[m] = 0;
					}
				}
			}
			ku1 << k << '\t' << b[k] << endl;
		}
		if (right == true)
		{
			i = nX - 1;
			k = NumberNode(xNet[i], yNet[j]);
			di[k] = 1;
			b[k] = AnaliticSolution(xNet[i], yNet[j]);
			for (m = ig[k]; m < ig[k + 1]; m++)
			{
				ggl[m] = 0;
			}
			for (l = 0; l < countOfNodes; l++)
			{
				for (m = ig[l]; m < ig[l + 1]; m++)
				{
					if (k == jg[m])
					{
						ggu[m] = 0;
					}
				}
			}
			ku1 << k << '\t' << b[k] << endl;
		}
	}
}

void Edge2(bool up, bool down, bool left, bool right)
{
	ofstream ku2("ku2.txt");
	int i, j, k1, k2, numbKE;
	double tetta1, tetta2, h = 0.01, dh;
	if (up == true)
	{
		j = nY - 1;
		for (i = 0; i < nX - 1; i++)
		{
			tetta1 = (AnaliticSolution(xNet[i], yNet[j] + h) - AnaliticSolution(xNet[i], yNet[j] - h)) / (2.0 * h);
			k1 = NumberNode(xNet[i], yNet[j]);
			tetta2 = (AnaliticSolution(xNet[i + 1], yNet[j] + h) - AnaliticSolution(xNet[i + 1], yNet[j] - h)) / (2.0 * h);
			k2 = NumberNode(xNet[i + 1], yNet[j]);
			dh = xNet[i + 1] - xNet[i];
			numbKE = FindNumberOfKe(k1, k2);
			b[k1] += dh * (2 * tetta1 + tetta2) / 6 * Lambda(numbKE);
			b[k2] += dh * (tetta1 + 2 * tetta2) / 6 * Lambda(numbKE);
			ku2 << setw(10) << k1 << setw(15) << tetta1 * Lambda(numbKE) << endl;
			ku2 << setw(10) << k2 << setw(15) << tetta2 * Lambda(numbKE) << endl;
		}
	}
	if (down == true)
	{
		j = 0;
		for (i = 0; i < nX - 1; i++)
		{
			tetta1 = -(AnaliticSolution(xNet[i], yNet[j] + h) - AnaliticSolution(xNet[i], yNet[j] - h)) / (2.0 * h);
			k1 = NumberNode(xNet[i], yNet[j]);
			dh = xNet[i + 1] - xNet[i];
			tetta2 = -(AnaliticSolution(xNet[i + 1], yNet[j] + h) - AnaliticSolution(xNet[i + 1], yNet[j] - h)) / (2.0 * h);
			k2 = NumberNode(xNet[i + 1], yNet[j]);
			numbKE = FindNumberOfKe(k1, k2);
			b[k1] += dh * (2 * tetta1 + tetta2) / 6 * Lambda(numbKE);
			b[k2] += dh * (tetta1 + 2 * tetta2) / 6 * Lambda(numbKE);
			ku2 << setw(10) << k1 << setw(15) << tetta1 * Lambda(numbKE) << endl;
			ku2 << setw(10) << k2 << setw(15) << tetta2 * Lambda(numbKE) << endl;
		}
	}

	if (left == true)
	{
		i = 0;
		for (j = 0; j < nY - 1; j++)
		{
			tetta1 = -(AnaliticSolution(xNet[i] + h, yNet[j]) - AnaliticSolution(xNet[i] - h, yNet[j])) / (2.0 * h);
			k1 = NumberNode(xNet[i], yNet[j]);
			dh = yNet[j + 1] - yNet[j];
			tetta2 = -(AnaliticSolution(xNet[i] + h, yNet[j + 1]) - AnaliticSolution(xNet[i] - h, yNet[j + 1])) / (2 * h);
			k2 = NumberNode(xNet[i], yNet[j + 1]);

			numbKE = FindNumberOfKe(k1, k2);
			b[k1] += dh * (2 * tetta1 + tetta2) / 6 * Lambda(numbKE);
			b[k2] += dh * (tetta1 + 2 * tetta2) / 6 * Lambda(numbKE);
			ku2 << setw(10) << k1 << setw(15) << tetta1 * Lambda(numbKE) << endl;
			ku2 << setw(10) << k2 << setw(15) << tetta2 * Lambda(numbKE) << endl;
		}
	}
	if (right == true)
	{
		i = nX - 1;
		for (j = 0; j < nY - 1; j++)
		{
			tetta1 = (AnaliticSolution(xNet[i] + h, yNet[j]) - AnaliticSolution(xNet[i] - h, yNet[j])) / (2.0 * h);
			k1 = NumberNode(xNet[i], yNet[j]);
			dh = yNet[j + 1] - yNet[j];
			tetta2 = (AnaliticSolution(xNet[i] + h, yNet[j + 1]) - AnaliticSolution(xNet[i] - h, yNet[j + 1])) / (2 * h);
			k2 = NumberNode(xNet[i], yNet[j + 1]);
			numbKE = FindNumberOfKe(k1, k2);
			b[k1] += dh * (2 * tetta1 + tetta2) / 6 * Lambda(numbKE);
			b[k2] += dh * (tetta1 + 2 * tetta2) / 6 * Lambda(numbKE);
			ku2 << setw(10) << k1 << setw(15) << tetta1 * Lambda(numbKE) << endl;
			ku2 << setw(10) << k2 << setw(15) << tetta2 * Lambda(numbKE) << endl;
		}
	}
}

void Edge3(bool up, bool down, bool left, bool right)
{
	ofstream ku3("ku3.txt");
	int i, j, k, k1, k2, numbKE, nextNode;
	double tetta1, tetta2, u, h = 0.01, dh, beta, ubeta1, ubeta2, leftNodeX, leftNodeY;
	if (up == true)
	{
		j = nY - 1;
		for (i = 0; i < nX - 1; i++)
		{
			leftNodeX = xNet[i];
			leftNodeY = yNet[j];
			k1 = NumberNode(xNet[i], yNet[j]);
			k2 = NumberNode(xNet[i + 1], yNet[j]);
			numbKE = FindNumberOfKe(k1, k2);
			tetta1 = (AnaliticSolution(leftNodeX, leftNodeY + h) - AnaliticSolution(leftNodeX, leftNodeY - h)) / (2.0 * h) * Lambda(numbKE);
			tetta2 = (AnaliticSolution(xNet[i + 1], yNet[j] + h) - AnaliticSolution(xNet[i + 1], yNet[j] - h)) / (2.0 * h) * Lambda(numbKE);
			ubeta1 = AnaliticSolution(leftNodeX, leftNodeY) - 1;
			ubeta2 = AnaliticSolution(xNet[i + 1], yNet[j]) - 1;
			beta = -tetta1;
			dh = xNet[i + 1] - leftNodeX;
			b[k1] += beta * dh * (2 * ubeta1 + ubeta2) / 6;
			b[k2] += beta * dh * (ubeta1 + 2 * ubeta2) / 6;
			di[k1] += beta * dh / 3;
			di[k2] += beta * dh / 3;
			for (k = ig[k2]; k < ig[k2 + 1]; k++)
			{
				if (jg[k] == k1)
				{
					ggl[k] += beta * dh / 6;
					break;
				}
			}
			ku3 << setw(15) << k1 << setw(15) << ubeta1 << setw(15) << k2 << setw(15) << ubeta2 << setw(15) << beta << endl;
		}
	}
	if (down == true)
	{
		j = 0;
		for (i = 0; i < nX - 1; i++)
		{
			leftNodeX = xNet[i];
			leftNodeY = yNet[j];
			k1 = NumberNode(xNet[i], yNet[j]);
			k2 = NumberNode(xNet[i + 1], yNet[j]);
			numbKE = FindNumberOfKe(k1, k2);
			tetta1 = -(AnaliticSolution(leftNodeX, leftNodeY + h) - AnaliticSolution(leftNodeX, leftNodeY - h)) / (2.0 * h) * Lambda(numbKE);
			tetta2 = -(AnaliticSolution(xNet[i + 1], yNet[j] + h) - AnaliticSolution(xNet[i + 1], yNet[j] - h)) / (2.0 * h) * Lambda(numbKE);
			ubeta1 = AnaliticSolution(leftNodeX, leftNodeY) - 1;
			ubeta2 = AnaliticSolution(xNet[i + 1], yNet[j]) - 1;
			beta = -tetta1;
			dh = xNet[i + 1] - leftNodeX;
			b[k1] += beta * dh * (2 * ubeta1 + ubeta2) / 6;
			b[k2] += beta * dh * (ubeta1 + 2 * ubeta2) / 6;
			di[k1] += beta * dh / 3;
			di[k2] += beta * dh / 3;
			for (k = ig[k2]; k < ig[k2 + 1]; k++)
			{
				if (jg[k] == k1)
				{
					ggl[k] += beta * dh / 6;
					break;
				}
			}
			ku3 << setw(15) << k1 << setw(15) << ubeta1 << setw(15) << k2 << setw(15) << ubeta2 << setw(15) << beta << endl;
		}
	}

	if (left == true)
	{
		i = 0;
		for (j = 0; j < nY - 1; j++)
		{
			leftNodeX = xNet[i];
			leftNodeY = yNet[j];
			k1 = NumberNode(xNet[i], yNet[j]);
			k2 = NumberNode(xNet[i], yNet[j + 1]);
			numbKE = FindNumberOfKe(k1, k2);
			tetta1 = -(AnaliticSolution(leftNodeX + h, leftNodeY) - AnaliticSolution(leftNodeX - h, leftNodeY)) / (2.0 * h) * Lambda(numbKE);
			tetta2 = -(AnaliticSolution(xNet[i] + h, yNet[j + 1]) - AnaliticSolution(xNet[i] - h, yNet[j + 1])) / (2.0 * h) * Lambda(numbKE);
			ubeta1 = AnaliticSolution(leftNodeX, leftNodeY) - 1;
			ubeta2 = AnaliticSolution(xNet[i], yNet[j + 1]) - 1;
			beta = -tetta1;
			dh = yNet[j + 1] - leftNodeY;
			b[k1] += beta * dh * (2 * ubeta1 + ubeta2) / 6;
			b[k2] += beta * dh * (ubeta1 + 2 * ubeta2) / 6;
			di[k1] += beta * dh / 3;
			di[k2] += beta * dh / 3;
			for (k = ig[k2]; k < ig[k2 + 1]; k++)
			{
				if (jg[k] == k1)
				{
					ggl[k] += beta * dh / 6;
					break;
				}
			}
			ku3 << setw(15) << k1 << setw(15) << ubeta1 << setw(15) << k2 << setw(15) << ubeta2 << setw(15) << beta << endl;
		}
	}
	if (right == true)
	{
		i = nX - 1;
		for (j = 0; j < nY - 1; j++)
		{
			leftNodeX = xNet[i];
			leftNodeY = yNet[j];
			k1 = NumberNode(xNet[i], yNet[j]);
			k2 = NumberNode(xNet[i], yNet[j + 1]);
			numbKE = FindNumberOfKe(k1, k2);
			tetta1 = (AnaliticSolution(leftNodeX + h, leftNodeY) - AnaliticSolution(leftNodeX - h, leftNodeY)) / (2.0 * h) * Lambda(numbKE);
			tetta2 = (AnaliticSolution(xNet[i] + h, yNet[j + 1]) - AnaliticSolution(xNet[i] - h, yNet[j + 1])) / (2.0 * h) * Lambda(numbKE);
			ubeta1 = AnaliticSolution(leftNodeX, leftNodeY) - 1;
			ubeta2 = AnaliticSolution(xNet[i], yNet[j + 1]) - 1;
			beta = -tetta1;
			dh = yNet[j + 1] - leftNodeY;
			b[k1] += beta * dh * (2 * ubeta1 + ubeta2) / 6;
			b[k2] += beta * dh * (ubeta1 + 2 * ubeta2) / 6;
			di[k1] += beta * dh / 3;
			di[k2] += beta * dh / 3;
			for (k = ig[k2]; k < ig[k2 + 1]; k++)
			{
				if (jg[k] == k1)
				{
					ggl[k] += beta * dh / 6;
					break;
				}
			}
			ku3 << setw(15) << k1 << setw(15) << ubeta1 << setw(15) << k2 << setw(15) << ubeta2 << setw(15) << beta << endl;
		}
	}
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

	//Edge2(1, 1, 1, 1);
	//Edge3(1, 1, 1, 1);
	//Edge1_sim(1, 1, 1, 1);
	//PrintPlotMatrix(true);
	for (i = 0; i < ig[countOfNodes]; i++)
	{
		ggu[i] = ggl[i];
	}
	Edge1_not_sim(1, 1, 1, 1);
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

double ScalarMult(double* v1, double* v2)
{
	int i;
	double result;
	result = 0;
	for (i = 0; i < countOfNodes; i++)
	{
		result += v1[i] * v2[i];
	}
	return result;
}

//	Умножение матрицы на вектор
void MultMatrixOnVector(double* in, double* out)
{
	int i, j;
	double* out1;
	out1 = new double[countOfNodes];
	for (i = 0; i < countOfNodes; i++)
	{
		out1[i] = di[i] * in[i];
		for (j = ig[i]; j < ig[i + 1]; j++)
		{
			out1[i] += ggl[j] * in[jg[j]];
			out1[jg[j]] += ggu[j] * in[i];
		}
	}
	for (i = 0; i < countOfNodes; i++)
		out[i] = out1[i];
	delete[] out1;
}

void runMSG()
{
	int maxiter = 10000, i;
	double alfa, alfachisl, alfaznam, beta, betachisl, betaznam, checkE, epsMSG = 1e-16;
	double* r = new double[countOfNodes];
	double* s = new double[countOfNodes];
	double* z = new double[countOfNodes];
	double* rout = new double[countOfNodes];
	for (i = 0; i < countOfNodes; i++)
	{
		s[i] = rout[i] = r[i] = q[i] = z[i] = 0;
	}
	MultMatrixOnVector(q, r);
	for (i = 0; i < countOfNodes; i++)
	{
		r[i] = b[i] - r[i];
		z[i] = r[i];
	}
	checkE = sqrt(ScalarMult(r, r) / ScalarMult(b, b));
	for (iter = 0; iter < maxiter && checkE >= epsMSG; iter++)
	{
		alfachisl = ScalarMult(r, r);
		MultMatrixOnVector(z, s);
		alfaznam = ScalarMult(s, z);
		alfa = alfachisl / alfaznam;
		for (i = 0; i < countOfNodes; i++)
		{
			q[i] = q[i] + alfa * z[i];
			rout[i] = r[i] - alfa * s[i];
		}
		betachisl = ScalarMult(rout, rout);
		betaznam = ScalarMult(r, r);
		beta = betachisl / betaznam;
		for (i = 0; i < countOfNodes; i++)
		{
			z[i] = rout[i] + beta * z[i];
			r[i] = rout[i];
		}
		checkE = sqrt(ScalarMult(r, r) / ScalarMult(b, b));
	}
}

void runLOS()
{
	int maxiter = 10000, i;
	double alfa, alfachisl, alfaznam, beta, betachisl, betaznam, checkE, epsMSG = 1e-16, startNeviazka;
	double* r = new double[countOfNodes];
	double* s = new double[countOfNodes];
	double* z = new double[countOfNodes];
	double* p = new double[countOfNodes];
	double* rout = new double[countOfNodes];
	for (i = 0; i < countOfNodes; i++)
	{
		s[i] = rout[i] = r[i] = q[i] = z[i] = p[i] = 0;
	}
	MultMatrixOnVector(q, r);
	for (i = 0; i < countOfNodes; i++)
	{
		r[i] = b[i] - r[i];
		z[i] = r[i];
	}
	MultMatrixOnVector(z, p);
	checkE = sqrt(ScalarMult(r, r) / ScalarMult(b, b));
	//startNeviazka = checkE = ScalarMult(r, r);
	for (iter = 0; iter < maxiter && checkE >= epsMSG; iter++)
	{
		alfachisl = ScalarMult(p, r);
		alfaznam = ScalarMult(p, p);
		alfa = alfachisl / alfaznam;
		for (i = 0; i < countOfNodes; i++)
		{
			q[i] = q[i] + alfa * z[i];
			r[i] = r[i] - alfa * p[i];
		}
		MultMatrixOnVector(r, rout);
		betachisl = ScalarMult(p, rout);
		betaznam = ScalarMult(p, p);
		beta = -betachisl / betaznam;
		for (i = 0; i < countOfNodes; i++)
		{
			z[i] = r[i] + beta * z[i];
			p[i] = rout[i] + beta * p[i];
		}
		checkE = sqrt(ScalarMult(r, r) / ScalarMult(b, b));
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
		outFile << setw(15) << q[i] << "\t" << setw(15) << analitSolvInNode << "\t" << setw(15) << abs(q[i] - analitSolvInNode) << endl;
		chislitPogreshnosti += (q[i] - analitSolvInNode) * (q[i] - analitSolvInNode);
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
	runLOS();

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
