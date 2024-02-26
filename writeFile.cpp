#include <fstream>
#include <iostream>
#include <vector>
#include <stdio.h>
#include <string>
#include "writeFile.h"
//#include "zapis_object.h"
//#include "vivod.h"

void write(int step, std:: vector<Node>& nodes, std:: vector<Cell>& cells, std:: vector <Part> parts) {
    
    std::string path = "res_" + std:: to_string(step) + ".vtu";
    std:: ofstream fout(path);
    if (!fout.is_open()) {
        std::cout << "Ошибка открытия\n";
    } 
    else 
    {
        fout << "<?xml version=\"1.0\"?>\n" 
             << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n" 
             << "<UnstructuredGrid>\n" 
             << "<Piece NumberOfPoints=\"" << nodes.size() << "\" " << "NumberOfCells=\"" << cells.size() << "\">\n" 
             << "<PointData Scalars=\"SpeedY\">\n" 
             << "<DataArray type=\"Float32\" Name=\"SpeedY\" format=\"ascii\">\n";
        for (const Node& node : nodes) 
        {
            fout << node.vz << " ";
        }
        fout << "\n</DataArray>\n"
             << "<DataArray type=\"Float32\" Name=\"SpeedX\" format=\"ascii\">\n";
        for (const Node& node : nodes) 
        {
            fout << node.vr << " ";
        }
        fout << "\n</DataArray>\n"
            << "</PointData>\n"
            << "<CellData Scalars=\"Part\">\n"
            << "<DataArray type=\"Int32\" Name=\"Part\" format=\"ascii\">\n";
        for (const Cell& cell : cells) 
        {
            fout <<cell.PID << " ";
        }
        fout << "\n</DataArray>\n"
             << "<DataArray type=\"Float32\" Name=\"Intensive_napr\" format=\"ascii\">\n";
        for (const Cell& cell : cells) 
        {
            fout << cell.intensive_sig << " ";
        }
        fout << "\n</DataArray>\n"
             << "</CellData>\n"
             << "<Points>\n"
             << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\"/>\n";
        for (const Node& node : nodes) 
        {
            fout <<node.r << " " << node.z << " " << 0.0 << "\n";
        }
        fout << "</Points>\n"
             << "<Cells>\n"
             << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
        for (const Cell& cell : cells)
        {
            fout << cell.NID[0] - 1 << " " << cell.NID[1] - 1 << " " << cell.NID[2] - 1 << " ";
        }
        fout << "\n</DataArray>\n"
             << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
        size_t i = 0;
        while (i <= 3 * cells.size())
        {
            i += 3;
            fout << i << " ";
        }
        fout << "\n</DataArray>\n"
             << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
        for (unsigned int i = 0; i < cells.size(); ++i)
        {
            fout << 5 << " ";
        }
        fout << "</DataArray>\n" << "</Cells>\n" << "</Piece>\n"
            << "</UnstructuredGrid>\n" << "</VTKFile>\n";
    }
}
