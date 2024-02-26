#include<string>
#include<vector>
#include<memory>
#ifndef READ_FILE
#define READ_FILE
//Данные, полученные из текстового файла
extern int part_vel; //Парт к которому приложена скорость
extern double init_vz; //Скорость начальная
extern double kr; //Число Куранта
extern double lamdal; //Линейная псевдовязкость
extern double lamdak; //Квадратичная псевдовязкость
extern double T_lim; //Время расчета
extern double T_zap; //Шаг записи
constexpr auto PI = 3.1415926535;
struct Cell
{
	int YID; //ID ячейки
	int PID; //ID парта, к которому принадлежит ячейка
	std::vector<int> NID;//ID узлов, принадлежащих ячейки
	int boundCheck = -1;
	int boundClass = -1;
	double p; //Давление
	double c; //Скорость звука
	double qkomb, qrr, qzz, q00, qrz; //Комбинированная псевдовязкость
	double err, ezz, e00, erz; //Скорость изменения деформаций
	double delta_err, delta_ezz, delta_e00, delta_erz; //Приращения деформаций
	double Drr, Dzz, D00, Drz; //Девиаторы напряжений
	double sigma_rr, sigma_zz, sigma_00, sigma_rz, intensive_sig; //Компоненты напряжения
	double E; //Энергия
	double perimeter; // Периметр
	double A, A12; //Площадь ячейки
	double r0; // Координата х центра ячейки
	double z0; // Координата z центра ячейки
	double M; //Масса ячейки
	double V, V12, V12tV; //Объем ячейки
	double ro0; //Начальная плотность ячейки
	double ro; //Плотность ячейки
	double roAV0;
};
struct Node
{
	int UID; //ID узла
	int checkPress = 0;
	int boundCheck = -1;
	double r, z; //Старые координаты узла
	double rn, zn; //Новые координаты узла
	double vr, vz; //Старая скорость узла
	double vrn, vzn; //Новая скорость узла
	std::vector<int> idCells; //ID ячеек, в которых числится узел
};
struct Part
{
	int id;
	int n_mat;
	std::vector<int> bondaryCells; //Граничные ячейки парта
};
class Mat
{
public:
	int id;
	double ro; //Плотность
	double Eu; //Модуль Юнга
	double kpuas; //Коэффициенет Пуассона
	double sigmay; //Предел текучести
	double G;
	double gamm[100];
	double aa[100], zz[100];
	double aa1[100], zz1[100];
};
class Me : public Mat
{
public:
	Me();
	void urs(int nmat, const Cell& cell, Cell& pseucell);
	void rachetcoeff(Cell& cell);
};
struct Pressure
{
	std::vector<double> time;
	std::vector<double> press;
};

std::pair<double, std::vector<double>> Perimeter(const std::vector<int>& nodeCell, const std::vector<Node>& nodes);
void InitGeometry(Cell& cell, std::vector<Part>& parts, std::vector<Mat>& mats, std::vector<Node>& nodes);
void SkipLine(std::string& s, std::istream& fin, int count);
void SkipValue(std::string& s, std::istream& fin, int count);
bool CheckBound(const std::vector<int>& nodeCell, std::vector<Node>& nodes, std::vector<Cell>& cells);
bool CheckAxis(const std::vector<int>& nodeCell, std::vector<Node>& nodes);
void BoundClass(Cell& cell, std::vector<Node>& nodes);
void ReadFile(std:: vector<Part>& Pr, std:: vector<Mat>& Mt, std:: vector<Cell>& Cells, std:: vector<Node>& Nodes, Pressure& Press);
void PostRead(std::vector<Part>& parts, std::vector<Mat>& mats, std::vector<Cell>& cells, std::vector<Node>& nodes);
#endif // READ_FILE