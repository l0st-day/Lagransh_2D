#include "readFile.h"
#ifndef EQUATIONS
#define EQUATIONS

extern double deltaT; //Шаг по времени
extern double deltaT12; //Половинный шаг по времени

double ApproxNode(double vel1, double vel2, double coord1, double coord2);
double ApproxCell(double stress1, double stress2, const Node& node1, const Node& node2);
void Geometry(const Cell& cell, Cell& pseucell, const std::vector<Node>& nodes);
void Kurant(const Cell& cell);
void FreeNodesSC(Node& node, const std::vector<Node>& nodes, std::vector<Cell>& cells);
void SpeedCoord(double time, std::vector<Node>& nodes, std::vector<Cell>& cells, std::vector<Part>& parts, Pressure& press);
void Defor(const Cell& cell, Cell& pseucell, std::vector<Node>& nodes, const double time);
void PseudoVis(const Cell& cell, Cell& pseucell);
void Stress(const Cell& cell, Cell& pseucell, std::vector<Node>& nodes, std::vector<Part>& parts, std::vector<Mat>& mats);
void Energy(Cell& cell, Cell& pseucell);

#endif // EQUATIONS

