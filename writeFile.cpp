// #include "windows.h"
// #include <stdlib.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <stdio.h>
#include <string>
#include "uravnenia_perep_cont_3.1.cpp"
//#include "zapis_object.h"
//#include "vivod.h"

void vivod(int step, vector<Uzel>& Uz, vector<Yacheika>& Ych, vector <Part> Pr) {
    
    string path = "res_" + to_string(step) + ".vtu";
    ofstream fout(path);
    if (!fout.is_open()) {
        cout << "Ошибка открытия" << endl;
    } else {
        fout << "<?xml version=\"1.0\"?>" << "\n";
        fout << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << "\n";
        fout << "<UnstructuredGrid>" << "\n";
        fout << "<Piece NumberOfPoints=\"" << Uz.size() << "\" " << "NumberOfCells=\"" << Ych.size() << "\">" << "\n";
        fout << "<PointData Scalars=\"SpeedY\">" << "\n";
        fout << "<DataArray type=\"Float32\" Name=\"SpeedY\" format=\"ascii\">" << "\n";
            for (unsigned int i = 0; i < Uz.size(); i++) {
                fout << Uz[i].vz << " ";
            }
        fout << "\n";
        fout << "</DataArray>" << "\n";
        fout << "<DataArray type=\"Float32\" Name=\"SpeedX\" format=\"ascii\">" << "\n";
            for (unsigned int i = 0; i < Uz.size(); i++) {
                fout << Uz[i].vr << " ";
            }
        fout << "\n";
        fout << "</DataArray>" << "\n";
        fout << "<DataArray type=\"Int32\" Name=\"Press_check\" format=\"ascii\">" << "\n";
            for (unsigned int i = 0; i < Uz.size(); i++) {
                fout << Uz[i].check_press << " ";
            }
        fout << "\n";
        fout << "</DataArray>" << "\n";
        fout << "<DataArray type=\"Int32\" Name=\"Gr_uz\" format=\"ascii\">" << "\n";
            for (unsigned int i = 0; i < Uz.size(); i++) {
                fout << Uz[i].Gr_uz << " ";
            }
        fout << "\n";
        fout << "</DataArray>" << "\n";
        fout << "<DataArray type=\"Int32\" Name=\"Ych_cont\" format=\"ascii\">" << "\n";
            for (unsigned int i = 0; i < Uz.size(); i++) {
                fout << Uz[i].Ych_cont << " ";
            }
        fout << "\n";
        fout << "</DataArray>" << "\n";
        fout << "</PointData>" << "\n";
        fout << "<CellData Scalars=\"Part\">" << "\n";
        fout << "<DataArray type=\"Int32\" Name=\"Part\" format=\"ascii\">" << "\n";
            for (unsigned int i = 0; i < Ych.size(); i++) {
                fout << Ych[i].PID << " ";
            }
        fout << "\n";
        fout << "</DataArray>" << "\n";
        fout << "<DataArray type=\"Float32\" Name=\"Intensive_napr\" format=\"ascii\">" << "\n";
            for (unsigned int i = 0; i < Ych.size(); i++) {
                fout << Ych[i].intensive_sig << " ";
            }
        fout << "\n";
        fout << "</DataArray>" << "\n";
        fout << "<DataArray type=\"Int32\" Name=\"Gran\" format=\"ascii\">" << "\n";
            for (unsigned int i = 0; i < Ych.size(); i++) {
                fout << Ych[i].Gr << " ";
            }
        fout << "\n";
        fout << "</DataArray>" << "\n";
        fout << "<DataArray type=\"Int32\" Name=\"Gran_klass\" format=\"ascii\">" << "\n";
            for (unsigned int i = 0; i < Ych.size(); i++) {
                fout << Ych[i].Gr_kl << " ";
            }
        fout << "\n";
        fout << "</DataArray>" << "\n";
        fout << "</CellData>" << "\n";
        fout << "<Points>" << "\n";
        fout << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\"/>" << "\n";
        for (unsigned int i = 0; i < Uz.size(); i++) {
            fout <<Uz[i].r << " " << Uz[i].z << " " << 0.0 << "\n";
        }
        // fout << "\n";
        // fout << "</DataArray>" << "\n";
        fout << "</Points>" << "\n";
        fout << "<Cells>" << "\n";
        fout << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">" << "\n";
            for (unsigned int i = 0; i < Ych.size(); i++)
            {
                fout << Ych[i].N1ID - 1 << " " << Ych[i].N2ID - 1 << " " << Ych[i].N3ID - 1 << " ";
            }
        fout << "\n";
        fout << "</DataArray>" << "\n";
        fout << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">" << "\n";
        unsigned int i = 0;
        while (i <= 3*Ych.size())
        {
            i += 3;
            fout << i << " ";
        }
        fout << "\n";
        fout << "</DataArray>" << "\n";
        fout << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">" << "\n";
            for (unsigned int i = 0; i < Ych.size(); i++)
            {
                fout << 5 << " ";
            }
        fout << "</DataArray>" << "\n";
        fout << "</Cells>" << "\n";
        fout << "</Piece>" << "\n";
        fout << "</UnstructuredGrid>" << "\n";
        fout << "</VTKFile>" << "\n";
    }
    fout.close();
}
