#include <cmath> 
#include "zapis_object.cpp"

//Треугольные ячейки =========================

double deltaT = 100.; //Шаг по времени
double deltaT12; //Половинный шаг по времени

//===================================
// Функциz апрокссимации узлов в центр ячейки
//===================================
double ApproxNode(double vel1, double vel2, double coord1, double coord2)
{
	return (vel2 + vel1) * (coord2 - coord1);
}
//Функции апрокссимации ячеек в узел
double ApproxCell(double stress1, double stress2, const Node& node1, const Node& node2) 
{
    return (stress1 * (node2.z - node1.z) - stress2 * (node2.r - node1.r));
}

int normal(std:: vector<Uzel>& Uz, int fUID, int aUID, int bUID, double& xf1, double& yf1, double& lff1)
{
	double bab = Uz[aUID].z / (1 + (Uz[bUID].z - Uz[aUID].z) * Uz[aUID].r / (Uz[aUID].z * Uz[bUID].r - Uz[bUID].z * Uz[aUID].r));
	double kab = bab * (Uz[bUID].z - Uz[aUID].z) / (Uz[aUID].z * Uz[bUID].r - Uz[bUID].z * Uz[aUID].r);
	double kff = -1.0 / kab;
	double bff = Uz[fUID].z - kff * Uz[fUID].r;
	xf1 = (bff - bab) / (kab - kff);
	yf1 = kff * xf1 + bff;
	lff1 = sqrt(pow(Uz[fUID].r - xf1, 2.0) + pow(Uz[fUID].z - yf1, 2.0));
	return 1;
};
//++++++++++++++++++++++++++++++++++++
int contact(int pos, std:: vector<Yacheika>& Ych, std:: vector<Uzel>& Uz, unsigned int i, int aUID, int bUID, double& lab, double& del,
	double& sinab, double& cosab, double& sigmad, int& id_ych_gr, double& chislitr, double& chislitz, double& omega, double& gama, double& beta, double& Gx, double& Gy, double& Rx, double& Ry, double& b, int& flag, Pressure& Press, double time)
{
	int eUID = 0;
	int id_ych = 0;
	double r1 = 0.0, r2 = 0.0, z1 = 0.0, z2 = 0.0;
	switch (Uz[i].kolvoY.size())
	{
	case 1:
	{
		id_ych = Uz[i].kolvoY[0] - 1;
		switch (pos)
		{
		case 1:
			if (Ych[id_ych].N2ID == Uz[i].UID)
			{
				eUID = Ych[id_ych].N1ID - 1;
				r1 = Uz[Ych[id_ych].N3ID - 1].r;
				r2 = Uz[Ych[id_ych].N1ID - 1].r;
				z1 = Uz[Ych[id_ych].N3ID - 1].z;
				z2 = Uz[Ych[id_ych].N1ID - 1].z;
			}
			else if (Ych[id_ych].N1ID == Uz[i].UID)
			{
				eUID = Ych[id_ych].N3ID - 1;
				r1 = Uz[Ych[id_ych].N2ID - 1].r;
				r2 = Uz[Ych[id_ych].N3ID - 1].r;
				z1 = Uz[Ych[id_ych].N2ID - 1].z;
				z2 = Uz[Ych[id_ych].N3ID - 1].z;
			}
			break;
		case -1:
			if (Ych[id_ych].N2ID == Uz[i].UID)
			{
				eUID = Ych[id_ych].N3ID - 1;
				r1 = Uz[Ych[id_ych].N3ID - 1].r;
				r2 = Uz[Ych[id_ych].N1ID - 1].r;
				z1 = Uz[Ych[id_ych].N3ID - 1].z;
				z2 = Uz[Ych[id_ych].N1ID - 1].z;
			}
			else if (Ych[id_ych].N1ID == Uz[i].UID)
			{
				eUID = Ych[id_ych].N2ID - 1;
				r1 = Uz[Ych[id_ych].N2ID - 1].r;
				r2 = Uz[Ych[id_ych].N3ID - 1].r;
				z1 = Uz[Ych[id_ych].N2ID - 1].z;
				z2 = Uz[Ych[id_ych].N3ID - 1].z;
			}
			break;
		default:
			break;
		}
		break;
	}
	default:
	{
		switch (pos)
		{
		case 1:
			for (unsigned int l = 0; l < Uz[i].kolvoY.size(); l++)
			{
				id_ych = Uz[i].kolvoY[l] - 1;
				if (Ych[id_ych].N2ID == Uz[i].UID && Uz[Ych[id_ych].N1ID - 1].Gr_uz == 1)
				{
					eUID = Ych[id_ych].N1ID - 1;
					r1 = Uz[Ych[id_ych].N3ID - 1].r;
					r2 = Uz[Ych[id_ych].N1ID - 1].r;
					z1 = Uz[Ych[id_ych].N3ID - 1].z;
					z2 = Uz[Ych[id_ych].N1ID - 1].z;
					break;
				}
				else if (Ych[id_ych].N3ID == Uz[i].UID && Uz[Ych[id_ych].N2ID - 1].Gr_uz == 1)
				{
					eUID = Ych[id_ych].N2ID - 1;
					r1 = Uz[Ych[id_ych].N1ID - 1].r;
					r2 = Uz[Ych[id_ych].N2ID - 1].r;
					z1 = Uz[Ych[id_ych].N1ID - 1].z;
					z2 = Uz[Ych[id_ych].N2ID - 1].z;
					break;
				}
			}
			break;
		case -1:
			for (unsigned int l = 0; l < Uz[i].kolvoY.size(); l++)
			{
				id_ych = Uz[i].kolvoY[l] - 1;
				if (Ych[id_ych].N1ID == Uz[i].UID && Uz[Ych[id_ych].N2ID - 1].Gr_uz == 1)
				{
					eUID = Ych[id_ych].N2ID - 1;
					r1 = Uz[Ych[id_ych].N2ID - 1].r;
					r2 = Uz[Ych[id_ych].N3ID - 1].r;
					z1 = Uz[Ych[id_ych].N2ID - 1].z;
					z2 = Uz[Ych[id_ych].N3ID - 1].z;
					break;
				}
				else if (Ych[id_ych].N2ID == Uz[i].UID && Uz[Ych[id_ych].N3ID - 1].Gr_uz == 1)
				{
					eUID = Ych[id_ych].N3ID - 1;
					r1 = Uz[Ych[id_ych].N3ID - 1].r;
					r2 = Uz[Ych[id_ych].N1ID - 1].r;
					z1 = Uz[Ych[id_ych].N3ID - 1].z;
					z2 = Uz[Ych[id_ych].N1ID - 1].z;
					break;
				}
				else if (Ych[id_ych].N3ID == Uz[i].UID && Uz[Ych[id_ych].N1ID - 1].Gr_uz == 1)
				{
					eUID = Ych[id_ych].N1ID - 1;
					r1 = Uz[Ych[id_ych].N1ID - 1].r;
					r2 = Uz[Ych[id_ych].N2ID - 1].r;
					z1 = Uz[Ych[id_ych].N1ID - 1].z;
					z2 = Uz[Ych[id_ych].N2ID - 1].z;
					break;
				}
			}
			break;

		default:
			break;
		}
		break;
	}
	}
	double xe1 = 0.0, ye1 = 0.0, lee1 = 0.0;
	normal(Uz, eUID, aUID, bUID, xe1, ye1, lee1);
	double lae1 = sqrt(pow(Uz[aUID].r - xe1, 2.0) + pow(Uz[aUID].z - ye1, 2.0));
	double lbe1 = sqrt(pow(Uz[bUID].r - xe1, 2.0) + pow(Uz[bUID].z - ye1, 2.0));
	// if (Uz[i].UID == 8025 && id_ych_gr == 7826 - 1)
	// {
	// 	cout << xe1 << " " << ye1 << " " << lee1 << " " << del << "\n";
	// 	cout << lae1 << " " << lbe1 << " " << lab - lbe1 - lae1 << "\n";
	// }
	if (fabs(lab - lae1 - lbe1) <= 1e-8)
	{
		if (lee1 <= del)
		{
			double lef = sqrt(pow(Uz[eUID].r - Uz[i].r, 2.0) + pow(Uz[eUID].z - Uz[i].z, 2.0));
			double sinef = 0.0, cosef = 0.0;
			chislitr = 0.0; chislitz = 0.0;
			switch (pos)
			{
			case 1:
				sinef = (Uz[eUID].z - Uz[i].z) / lef;
				cosef = (Uz[eUID].r - Uz[i].r) / lef;
				break;
			case -1:
				sinef = (Uz[i].z - Uz[eUID].z) / lef;
				cosef = (Uz[i].r - Uz[eUID].r) / lef;
				break;
			default:
				break;
			}
			sinab = (Uz[bUID].z - Uz[aUID].z) / lab;
			cosab = (Uz[bUID].r - Uz[aUID].r) / lab;
			sigmad = (Ych[id_ych_gr].sigma_rr * pow(sinef, 2.0) + Ych[id_ych_gr].sigma_zz * pow(cosef, 2.0) -
				2.0 * Ych[id_ych_gr].sigma_rz * sinef * cosef);
			b += lef * Ych[id_ych_gr].roAV0 / lab;
			chislitr = (aproxSig_r(Ych[id_ych].sigma_rr, Ych[id_ych].sigma_rz, r1, r2, z1, z2) - sigmad * sinef * lef) / (2.0 * omega);
			chislitz = (aproxSig_z(Ych[id_ych].sigma_zz, Ych[id_ych].sigma_rz, r1, r2, z1, z2) + sigmad * cosef * lef) / (2.0 * omega);
			// if (Uz[i].UID == 8025 && id_ych_gr == 7826 - 1)
			// {
			// 	cout << "Добавки от контакта" << "\n";
			// 	cout << lef << " " << sinef << " " << cosef << "\n";
			// 	cout << sinab << " " << cosab << "\n";
			// 	cout << sigmad << " " << b << " " << chislitr << " " << chislitz << "\n";
			// }
			flag += 1;
		}
	}
	else
	{
		switch (pos)
		{
		case 1:
		{
			double xb1 = 0.0, yb1 = 0.0, lbb1 = 0.0;
			normal(Uz, bUID, i, eUID, xb1, yb1, lbb1);
			if (lbb1 <= del)
			{
				double lef = sqrt(pow(Uz[eUID].r - Uz[i].r, 2.0) + pow(Uz[eUID].z - Uz[i].z, 2.0));
				double lfb = sqrt(pow(Uz[bUID].r - Uz[i].r, 2.0) + pow(Uz[bUID].z - Uz[i].z, 2.0));
				double sinef = (Uz[eUID].z - Uz[i].z) / lef;
				double cosef = (Uz[eUID].r - Uz[i].r) / lef;
				chislitr = 0.0; chislitz = 0.0;
				sinab = (Uz[bUID].z - Uz[aUID].z) / lab;
				cosab = (Uz[bUID].r - Uz[aUID].r) / lab;
				sigmad = lfb * (Ych[id_ych_gr].sigma_rr * pow(sinef, 2.0) + Ych[id_ych_gr].sigma_zz * pow(cosef, 2.0) -
					2.0 * Ych[id_ych_gr].sigma_rz * sinef * cosef) / lef;
				b += lfb * Ych[id_ych_gr].roAV0 / lab;
				for (unsigned int m = 0; m < Uz[bUID].kolvoY.size(); m++)
				{
					if ((Ych[Uz[bUID].kolvoY[m] - 1].N1ID - 1 == bUID &&
						Uz[Ych[Uz[bUID].kolvoY[m] - 1].N2ID - 1].Gr_uz == 1) ||
						(Ych[Uz[bUID].kolvoY[m] - 1].N2ID - 1 == bUID &&
							Uz[Ych[Uz[bUID].kolvoY[m] - 1].N3ID - 1].Gr_uz == 1))
					{
						double dUID = 0;
						if (Ych[Uz[bUID].kolvoY[m] - 1].N1ID - 1 == bUID)
						{
							dUID = Ych[Uz[bUID].kolvoY[m] - 1].N2ID - 1;
						}
						else if (Ych[Uz[bUID].kolvoY[m] - 1].N2ID - 1 == bUID)
						{
							dUID = Ych[Uz[bUID].kolvoY[m] - 1].N3ID - 1;
						}
						xe1 = 0.0, ye1 = 0.0, lee1 = 0.0;
						normal(Uz, eUID, bUID, dUID, xe1, ye1, lee1);
						if (lee1 <= del)
						{
							double lbe = sqrt(pow(Uz[eUID].r - Uz[bUID].r, 2.0) + pow(Uz[eUID].z - Uz[bUID].z, 2.0));
							double lbd = sqrt(pow(Uz[dUID].r - Uz[bUID].r, 2.0) + pow(Uz[dUID].z - Uz[bUID].z, 2.0));
							sigmad += lbe * (Ych[Uz[bUID].kolvoY[m] - 1].sigma_rr * pow(sinef, 2.0) + Ych[Uz[bUID].kolvoY[m] - 1].sigma_zz * pow(cosef, 2.0) -
								2.0 * Ych[Uz[bUID].kolvoY[m] - 1].sigma_rz * sinef * cosef) / lef;
							b += lbe * Ych[Uz[bUID].kolvoY[m] - 1].roAV0 / lbd;
						}
					}
				}
				chislitr = (aproxSig_r(Ych[id_ych].sigma_rr, Ych[id_ych].sigma_rz, r1, r2, z1, z2) - sigmad * (Uz[eUID].z - Uz[i].z)) / (2.0 * omega);
				chislitz = (aproxSig_z(Ych[id_ych].sigma_zz, Ych[id_ych].sigma_rz, r1, r2, z1, z2) + sigmad * (Uz[eUID].r - Uz[i].r)) / (2.0 * omega);
				flag += 1;
			}
			break;
		}
		case -1:
		{
			double xa1 = 0.0, ya1 = 0.0, laa1 = 0.0;
			normal(Uz, aUID, eUID, i, xa1, ya1, laa1);
			if (laa1 <= del)
			{
				double lef = sqrt(pow(Uz[eUID].r - Uz[i].r, 2.0) + pow(Uz[eUID].z - Uz[i].z, 2.0));
				double lfa = sqrt(pow(Uz[aUID].r - Uz[i].r, 2.0) + pow(Uz[aUID].z - Uz[i].z, 2.0));
				double sinef = (Uz[i].z - Uz[eUID].z) / lef;
				double cosef = (Uz[i].r - Uz[eUID].r) / lef;
				chislitr = 0.0; chislitz = 0.0;
				sinab = (Uz[bUID].z - Uz[aUID].z) / lab;
				cosab = (Uz[bUID].r - Uz[aUID].r) / lab;
				sigmad = lfa * (Ych[id_ych_gr].sigma_rr * pow(sinef, 2.0) + Ych[id_ych_gr].sigma_zz * pow(cosef, 2.0) -
					2.0 * Ych[id_ych_gr].sigma_rz * sinef * cosef) / lef;
				b += lfa * Ych[id_ych_gr].roAV0 / lab;
				for (unsigned int m = 0; m < Uz[aUID].kolvoY.size(); m++)
				{
					if ((Ych[Uz[aUID].kolvoY[m] - 1].N2ID - 1 == aUID &&
						Uz[Ych[Uz[aUID].kolvoY[m] - 1].N1ID - 1].Gr_uz == 1) ||
						(Ych[Uz[aUID].kolvoY[m] - 1].N3ID - 1 == aUID &&
							Uz[Ych[Uz[aUID].kolvoY[m] - 1].N2ID - 1].Gr_uz == 1))
					{
						double dUID = 0;
						if (Ych[Uz[aUID].kolvoY[m] - 1].N2ID - 1 == aUID)
						{
							dUID = Ych[Uz[aUID].kolvoY[m] - 1].N1ID - 1;
						}
						else if (Ych[Uz[aUID].kolvoY[m] - 1].N3ID - 1 == aUID)
						{
							dUID = Ych[Uz[aUID].kolvoY[m] - 1].N2ID - 1;
						}
						xe1 = 0.0, ye1 = 0.0, lee1 = 0.0;
						normal(Uz, eUID, dUID, aUID, xe1, ye1, lee1);
						if (lee1 <= del)
						{
							double lae = sqrt(pow(Uz[eUID].r - Uz[aUID].r, 2.0) + pow(Uz[eUID].z - Uz[aUID].z, 2.0));
							double lad = sqrt(pow(Uz[dUID].r - Uz[aUID].r, 2.0) + pow(Uz[dUID].z - Uz[aUID].z, 2.0));
							sigmad += lae * (Ych[Uz[aUID].kolvoY[m] - 1].sigma_rr * pow(sinef, 2.0) + Ych[Uz[aUID].kolvoY[m] - 1].sigma_zz * pow(cosef, 2.0) -
								2.0 * Ych[Uz[aUID].kolvoY[m] - 1].sigma_rz * sinef * cosef) / lef;
							b += lae * Ych[Uz[bUID].kolvoY[m] - 1].roAV0 / lad;
						}
					}
				}
				chislitr = (aproxSig_r(Ych[id_ych].sigma_rr, Ych[id_ych].sigma_rz, r1, r2, z1, z2) - sigmad * (Uz[i].z - Uz[eUID].z)) / (2.0 * omega);
				chislitz = (aproxSig_z(Ych[id_ych].sigma_zz, Ych[id_ych].sigma_rz, r1, r2, z1, z2) + sigmad * (Uz[i].r - Uz[eUID].r)) / (2.0 * omega);
				flag += 1;
			}
			break;
		}
		default:
			break;
		}
	}
	if (flag != 0)
	{
		for (unsigned int l = 0; l < Uz[i].kolvoY.size(); l++)
		{
			id_ych = Uz[i].kolvoY[l] - 1;
			switch (Ych[id_ych].Gr)
			{
			case -1:
				if (Uz[i].UID == Uz[Ych[id_ych].N1ID - 1].UID)
				{
					r1 = Uz[Ych[id_ych].N2ID - 1].r;
					r2 = Uz[Ych[id_ych].N3ID - 1].r;
					z1 = Uz[Ych[id_ych].N2ID - 1].z;
					z2 = Uz[Ych[id_ych].N3ID - 1].z;
				}
				else if (Uz[i].UID == Uz[Ych[id_ych].N2ID - 1].UID)
				{
					r1 = Uz[Ych[id_ych].N3ID - 1].r;
					r2 = Uz[Ych[id_ych].N1ID - 1].r;
					z1 = Uz[Ych[id_ych].N3ID - 1].z;
					z2 = Uz[Ych[id_ych].N1ID - 1].z;
				}
				else if (Uz[i].UID == Uz[Ych[id_ych].N3ID - 1].UID)
				{
					r1 = Uz[Ych[id_ych].N1ID - 1].r;
					r2 = Uz[Ych[id_ych].N2ID - 1].r;
					z1 = Uz[Ych[id_ych].N1ID - 1].z;
					z2 = Uz[Ych[id_ych].N2ID - 1].z;
				}
				chislitr += aproxSig_r(Ych[id_ych].sigma_rr, Ych[id_ych].sigma_rz, r1, r2, z1, z2) / (2.0 * omega);
				chislitz += aproxSig_z(Ych[id_ych].sigma_zz, Ych[id_ych].sigma_rz, r1, r2, z1, z2) / (2.0 * omega);
				break;
			default:
				break;
			}
			if (Uz[i].check_press > 0 && Ych[id_ych].Gr == 1)
				{
					/*if (Uz[i].UID == 7935)
					{
						cout << "PRESS\n";
						cout << Ych[id_ych].YID << "\n";
					}*/

					double r1_p = 0.0, r2_p = 0.0, z1_p = 0.0, z2_p = 0.0;
					if ((Ych[id_ych].N1ID == Uz[i].UID && Uz[Ych[id_ych].N2ID - 1].check_press > 0) ||
						(Ych[id_ych].N2ID == Uz[i].UID && Uz[Ych[id_ych].N1ID - 1].check_press > 0))
					{
						r1_p = Uz[Ych[id_ych].N1ID - 1].r;
						r2_p = Uz[Ych[id_ych].N2ID - 1].r;
						z1_p = Uz[Ych[id_ych].N1ID - 1].z;
						z2_p = Uz[Ych[id_ych].N2ID - 1].z;
						for (unsigned int k = 0; k < Press.time.size() - 1; k++)
						{
							if (time >= Press.time[k] && time <= Press.time[k + 1])
							{
								double press = Press.press[k] +
									((Press.press[k + 1] - Press.press[k]) * (time - Press.time[k])) /
									(Press.time[k + 1] - Press.time[k]);
								chislitz += press * (r2_p - r1_p) / sqrt(pow(r2_p - r1_p, 2) + pow(z2_p - z1_p, 2));
								chislitr -= press * (z2_p - z1_p) / sqrt(pow(r2_p - r1_p, 2) + pow(z2_p - z1_p, 2));
								/*if (Uz[i].UID == 7935)
								{
									cout << press << "\n";
									cout << chislitz << "\n";
									cout << chislitr << "\n";
								}*/
							}
						}
					}
					if ((Ych[id_ych].N3ID == Uz[i].UID && Uz[Ych[id_ych].N2ID - 1].check_press > 0) ||
						(Ych[id_ych].N2ID == Uz[i].UID && Uz[Ych[id_ych].N3ID - 1].check_press > 0))
					{
						r1_p = Uz[Ych[id_ych].N2ID - 1].r;
						r2_p = Uz[Ych[id_ych].N3ID - 1].r;
						z1_p = Uz[Ych[id_ych].N2ID - 1].z;
						z2_p = Uz[Ych[id_ych].N3ID - 1].z;
						for (unsigned int k = 0; k < Press.time.size() - 1; k++)
						{
							if (time >= Press.time[k] && time <= Press.time[k + 1])
							{
								double press = Press.press[k] +
									((Press.press[k + 1] - Press.press[k]) * (time - Press.time[k])) /
									(Press.time[k + 1] - Press.time[k]);
								chislitz += press * (r2_p - r1_p) / 2;
								chislitr -= press * (z2_p - z1_p) / 2;
								/*if (Uz[i].UID == 7935)
								{
									cout << press << "\n";
									cout << chislitz << "\n";
									cout << chislitr << "\n";
								}*/
							}
						}
					}
					
					/*if (Uz[i].UID == 7935)
					{
						cout << r1_p << " " << r2_p << " " << z1_p << " " << z2_p << "\n";
					}*/
				}
		}
		// if (Uz[i].UID == 8025)
		// {	
		// 	double g = (chislitr * cosab + chislitz * sinab) * cosab + 0.5 * (gama * cosab + beta * sinab) * cosab;
		// 	double r = (-chislitr * sinab + chislitz * cosab) * cosab + 0.5 * (-gama * sinab + beta * cosab) * cosab;
		// 	cout << id_ych_gr + 1 << "\n";
		// 	cout << g << " " << r << "\n";
		// }
		Gx += (chislitr * cosab + chislitz * sinab) * cosab + 0.5 * (gama * cosab + beta * sinab) * cosab;
		Gy += (chislitr * cosab + chislitz * sinab) * sinab + 0.5 * (gama * cosab + beta * sinab) * sinab;
		Rx += (chislitr * sinab + chislitz * cosab) * sinab + 0.5 * (gama * sinab + beta * cosab) * sinab;
		Ry += (chislitr * sinab + chislitz * cosab) * cosab + 0.5 * (gama * sinab + beta * cosab) * cosab;
	}
	return 1;
};
//+++++++++++++++++++++++



//===================================
//Расчет геометрии 
//===================================
void Geometry(const Cell& cell, Cell& pseucell, const std:: vector<Node>& nodes) 
{
	pseucell.perimeter = Perimeter(cell.NID, nodes).first;
	std:: vector<double> lengsid = Perimeter(cell.NID, nodes).second;
	double per = 0.5 * cell.perimeter;
	pseucell.A = sqrt(per * (per - lengsid[0]) * (per - lengsid[1]) * (per - lengsid[2]));
	for (size_t i = 0; i != cell.NID.size(); ++i)
	{
		pseucell.r0 += nodes[cell.NID[i] - 1].rn;
		pseucell.z0 += nodes[cell.NID[i] - 1].zn;
	}
	pseucell.r0 /= 3; pseucell.z0 /= 3;
	pseucell.V = 2 * PI * cell.ro0 * pseucell.A * pseucell.r0 / cell.M;
	pseucell.ro = cell.ro0 / pseucell.V;
}
//===================================
//УРС
//===================================
Me::Me()
{
	aa[5] = 3.02e+10; zz[5] = 4.8; gamm[5] = 0; //����
	aa[8] = 1.97e+10; zz[8] = 4.2; gamm[8] = 0; //��������
	aa[11] = 2.15e+10; zz[11] = 5.5; gamm[11] = 0; //�����
	aa[12] = 1.8797e9; zz[12] = 4.668; gamm[12] = 0; //����
	aa[13] = 3.595e8; zz[13] = 1.736; gamm[13] = 0; //���������� 320
	aa[14] = 0.86e+10; zz[14] = 5.3; gamm[14] = 0; //������
	aa[15] = 4.7353e6; zz[15] = 3.46; gamm[15] = 0; //���������� 160
	aa[17] = 1.4506e+11; zz[17] = 1.987; gamm[17] = 0; //������
	aa[19] = 1.955e+10; zz[19] = 3.648; gamm[19] = 0; //LiF
	//aa[19] = 1.472e+9; zz[19] = 8.76; gamm[19] = 0; //������ ��������� ��
	//aa[19] = 1.97e+10; zz[19] = 4.2; gamm[19] = 0;
	aa[21] = 5.71e+10; zz[21] = 4.7; gamm[21] = 0; // ���-90
	aa[22] = 2.3193e8; zz[22] = 3.539; gamm[22] = 0; //���������� 626
	aa[24] = 2.15e+10; zz[24] = 5.5; gamm[24] = 2.6; //�����
};
void Me::rachetcoeff(Cell& cell)
{
	for (int i = 5; i <= 99; i++)
	{
		zz1[i] = zz[i] - 1.0;
		aa1[i] = aa[i] * zz[i] / cell.ro0;
	}
};
void Me::urs(int nmat, const Cell& cell, Cell& pseucell)
{
	double a1 = aa[nmat];
	double z1 = zz[nmat];
	double az1 = aa1[nmat];
	double z11 = zz1[nmat];
	pseucell.p = a1 * (pow((pseucell.ro / cell.ro0), z1) - 1);

	double as2 = az1 * pow((pseucell.ro / cell.ro0), z11);
	if (as2 > 0)
	{
		pseucell.c = sqrt(as2);
	}
	else
	{
		pseucell.c = 100.;
	}
};
//===================================
//Расчет шага по времени
//===================================
void Kurant(const Cell& cell) 
{
    double deltat =  kr * 2.0 * cell.A / (cell.perimeter * cell.c);    
    // double f = 1 / (kr * kr);
    // double beta = 2.0 * f * (1.0 - f) / (1.0 - 2.0 * f);
    // double deltsht = (cell.perimetr / (3.0 * cell.c)) * sqrt(2.0 / (1.0 + beta + sqrt(1.0 + pow(beta, 2.0))));
    // if (deltsht < deltat)
    // {
    //     deltat = deltsht;
    // }
    if(deltat < deltaT) 
	{
        deltaT12 = (deltaT + deltat)/2.0;
        deltaT = deltat;
    }
};


void FreeNodesSC(Node& node, const std:: vector<Node>& nodes)
{
	double approx_r = 0.0,
			approx_z = 0.0,
			omega = 0.0,
			beta = 0.0,
			gama = 0.0;
	for(const auto cell : node.idCells)
	{
		if (cell -> NID[0] == node.UID)
		{
			approx_r += ApproxCell(cell -> sigma_rr, cell -> sigma_rz, nodes[cell -> NID[1] - 1], nodes[cell -> NID[2] - 1]);
			approx_z += ApproxCell(cell -> sigma_rz, cell -> sigma_zz, nodes[cell -> NID[1] - 1], nodes[cell -> NID[2] - 1]);
		}
		else if (cell -> NID[1] == node.UID)
		{
			approx_r += ApproxCell(cell -> sigma_rr, cell -> sigma_rz, nodes[cell -> NID[2] - 1], nodes[cell -> NID[1] - 1]);
			approx_z += ApproxCell(cell -> sigma_rz, cell -> sigma_zz, nodes[cell -> NID[2] - 1], nodes[cell -> NID[1] - 1]);
		}
		else
		{
			approx_r += ApproxCell(cell -> sigma_rr, cell -> sigma_rz, nodes[cell -> NID[1] - 1], nodes[cell -> NID[0] - 1]);
			approx_z += ApproxCell(cell -> sigma_rz, cell -> sigma_zz, nodes[cell -> NID[1] - 1], nodes[cell -> NID[0] - 1]);
		}
		omega += cell -> ro * cell -> A;
		beta += (cell -> sigma_rr - cell -> sigma_00) * cell -> A / cell -> M;
		gama += cell -> sigma_rz * cell -> A / cell -> M;
	}
	omega /= 3;
	beta /= node.idCells.size();
	gama /= node.idCells.size();
	node.vrn = node.vr + (0.5 * approx_r / omega + beta) * deltaT;
	node.vzn = node.vz + (0.5 * approx_z / omega + gama) * deltaT;
	node.rn = node.r + node.vrn * deltaT12;
	node.zn = node.z + node.vzn * deltaT12;
}


//===================================
//Уравнение движения
//===================================
void SpeedCoord(double time, std:: vector<Node>& nodes, std:: vector<Cell>& cells, std:: vector<Part> parts, std:: vector<Pressure> press)
{

	for (Node& node : nodes)
	{
		if (node.UID == -1)
		{
			node.r = 0.0;
			node.z = 0.0;
			node.vr = 0.0;
			node.vz = 0.0;
			continue;
		}
		
		switch (node.boundCheck)
		{
		case -1:
		{
			FreeNodesSC(node, nodes);
			break;
		}
		case 0:
		{
			FreeNodesSC(node, nodes);
			node.rn = 0.0;
			node.vrn = 0.0;
			break;
		}
		case 1:
		{
			int flad = 0, flag = 0;
			double b = 0.0;
			double z = 0.0;
			double Gx = 0.0;
			double Gy = 0.0;
			double Rx = 0.0;
			double Ry = 0.0;
			double l = 0.0;
			double m = 0.0;
			double deltaN = 0.0;
			double znamen_m = 0.0;
			for (unsigned int j = 0; j < Uz[i].kolvoY.size(); j++)
			{
				int id_ych = Uz[i].kolvoY[j] - 1;
				double r1 = 0.0, r2 = 0.0, z1 = 0.0, z2 = 0.0;
				if (Ych[id_ych].N1ID == Uz[i].UID)
				{
					r1 = Uz[Ych[id_ych].N2ID - 1].r;
					r2 = Uz[Ych[id_ych].N3ID - 1].r;
					z1 = Uz[Ych[id_ych].N2ID - 1].z;
					z2 = Uz[Ych[id_ych].N3ID - 1].z;
				}
				else if (Ych[id_ych].N2ID == Uz[i].UID)
				{
					r1 = Uz[Ych[id_ych].N3ID - 1].r;
					r2 = Uz[Ych[id_ych].N1ID - 1].r;
					z1 = Uz[Ych[id_ych].N3ID - 1].z;
					z2 = Uz[Ych[id_ych].N1ID - 1].z;
				}
				else if (Ych[id_ych].N3ID == Uz[i].UID)
				{
					r1 = Uz[Ych[id_ych].N1ID - 1].r;
					r2 = Uz[Ych[id_ych].N2ID - 1].r;
					z1 = Uz[Ych[id_ych].N1ID - 1].z;
					z2 = Uz[Ych[id_ych].N2ID - 1].z;
				};
				chislitr += aproxSig_r(Ych[id_ych].sigma_rr, Ych[id_ych].sigma_rz, r1, r2, z1, z2);
				chislitz += aproxSig_z(Ych[id_ych].sigma_zz, Ych[id_ych].sigma_rz, r1, r2, z1, z2);
				omega += Ych[id_ych].ro0 * Ych[id_ych].A / Ych[id_ych].V;
				beta += (Ych[id_ych].sigma_rr - Ych[id_ych].sigma_00) *
					(Ych[id_ych].A / Ych[id_ych].M);
				gama += Ych[id_ych].sigma_rz * (Ych[id_ych].A / Ych[id_ych].M);
				znamen_m += Ych[id_ych].roAV0;
				// if (Uz[i].UID == 6743)
				// {
				// 	cout << chislitr << " " << chislitz << " " << omega << " " << beta << " " << gama << " " << znamen_m << "\n";
				// }
				if (Uz[i].check_press > 0 && Ych[id_ych].Gr == 1)
				{
					/*if (Uz[i].UID == 7935)
					{
						cout << "PRESS\n";
						cout << Ych[id_ych].YID << "\n";
					}*/

					double r1_p = 0.0, r2_p = 0.0, z1_p = 0.0, z2_p = 0.0;
					if ((Ych[id_ych].N1ID == Uz[i].UID && Uz[Ych[id_ych].N2ID - 1].check_press > 0) ||
						(Ych[id_ych].N2ID == Uz[i].UID && Uz[Ych[id_ych].N1ID - 1].check_press > 0))
					{
						r1_p = Uz[Ych[id_ych].N1ID - 1].r;
						r2_p = Uz[Ych[id_ych].N2ID - 1].r;
						z1_p = Uz[Ych[id_ych].N1ID - 1].z;
						z2_p = Uz[Ych[id_ych].N2ID - 1].z;
						for (unsigned int k = 0; k < Press.time.size() - 1; k++)
						{
							if (time >= Press.time[k] && time <= Press.time[k + 1])
							{
								double press = Press.press[k] +
									((Press.press[k + 1] - Press.press[k]) * (time - Press.time[k])) /
									(Press.time[k + 1] - Press.time[k]);
								chislitz += press * (r2_p - r1_p) / 2;
								chislitr -= press * (z2_p - z1_p) / 2;
								/*if (Uz[i].UID == 7935)
								{
									cout << press << "\n";
									cout << chislitz << "\n";
									cout << chislitr << "\n";
								}*/
							}
						}
					}
					if ((Ych[id_ych].N3ID == Uz[i].UID && Uz[Ych[id_ych].N2ID - 1].check_press > 0) ||
						(Ych[id_ych].N2ID == Uz[i].UID && Uz[Ych[id_ych].N3ID - 1].check_press > 0))
					{
						r1_p = Uz[Ych[id_ych].N2ID - 1].r;
						r2_p = Uz[Ych[id_ych].N3ID - 1].r;
						z1_p = Uz[Ych[id_ych].N2ID - 1].z;
						z2_p = Uz[Ych[id_ych].N3ID - 1].z;
						for (unsigned int k = 0; k < Press.time.size() - 1; k++)
						{
							if (time >= Press.time[k] && time <= Press.time[k + 1])
							{
								double press = Press.press[k] +
									((Press.press[k + 1] - Press.press[k]) * (time - Press.time[k])) /
									(Press.time[k + 1] - Press.time[k]);
								chislitz += press * (r2_p - r1_p) / 2;
								chislitr -= press * (z2_p - z1_p) / 2;
								/*if (Uz[i].UID == 7935)
								{
									cout << press << "\n";
									cout << chislitz << "\n";
									cout << chislitr << "\n";
								}*/
							}
						}
					}
					
					/*if (Uz[i].UID == 7935)
					{
						cout << r1_p << " " << r2_p << " " << z1_p << " " << z2_p << "\n";
					}*/
				}
			}
			omega /= Uz[i].kolvoY.size();
			beta /= Uz[i].kolvoY.size();
			gama /= Uz[i].kolvoY.size();
			for (unsigned int j = 0; j < Pr.size(); j++)
			{
				if (Ych[Uz[i].kolvoY[0] - 1].PID != Pr[j].id)
				{
					for (unsigned k = 0; k < Pr[j].id_yach_gran.size(); k++)
					{
						int id_ych_gr = Pr[j].id_yach_gran[k] - 1;
						int aUID = 0, bUID = 0, cUID = 0;
						double sinab = 0.0, cosab = 0.0;
						double sigmad = 0.0;
						switch (Ych[id_ych_gr].Gr_kl)
						{
						case 0:
						{
							if (Uz[Ych[id_ych_gr].N1ID - 1].Gr_uz == 1 &&
								Uz[Ych[id_ych_gr].N2ID - 1].Gr_uz == 1)
							{
								aUID = Ych[id_ych_gr].N1ID - 1;
								bUID = Ych[id_ych_gr].N2ID - 1;
								cUID = Ych[id_ych_gr].N3ID - 1;
							}
							else if (Uz[Ych[id_ych_gr].N2ID - 1].Gr_uz == 1 &&
								Uz[Ych[id_ych_gr].N3ID - 1].Gr_uz == 1)
							{
								aUID = Ych[id_ych_gr].N2ID - 1;
								bUID = Ych[id_ych_gr].N3ID - 1;
								cUID = Ych[id_ych_gr].N1ID - 1;
							}
							double xf1 = 0.0, yf1 = 0.0, lff1 = 0.0;
							normal(Uz, i, aUID, bUID, xf1, yf1, lff1);
							double lab = sqrt(pow(Uz[aUID].r - Uz[bUID].r, 2.0) + pow(Uz[aUID].z - Uz[bUID].z, 2.0));
							double laf1 = sqrt(pow(Uz[aUID].r - xf1, 2.0) + pow(Uz[aUID].z - yf1, 2.0));
							double lbf1 = sqrt(pow(Uz[bUID].r - xf1, 2.0) + pow(Uz[bUID].z - yf1, 2.0));
							// if (Uz[i].UID == 60 && id_ych_gr == 3842 - 1)
							// {
							// 	cout << "Контакт найден" << "\n";
							// 	cout << aUID << " " << bUID << " " << cUID << "\n";
							// 	cout << lff1 << " " << lab << " " << 1e-4 / lab << "\n";
							// 	cout << lab - laf1 - lbf1 << "\n";
							// }
							// if (Uz[i].UID == 4861 && id_ych_gr == 13023)
							// {
							// 	cout << aUID + 1 << " " << bUID + 1 << "\n";
							// 	cout << lff1 << " " << 1e-14 / lab << " " << lab - laf1 - lbf1 << "\n";
							// }
							if (lff1 <= 1e-8 / lab && fabs(lab - laf1 - lbf1) <= 1e-6)
							{
								// if (Uz[i].UID == 60)
								// {
								// 	cout << "Контакт с ячейкой: " << id_ych_gr + 1 << "\n";
								// }
								double A = Uz[bUID].z - Uz[aUID].z;
								double B = -(Uz[bUID].r - Uz[aUID].r);
								if ((Uz[cUID].r - Uz[aUID].r) * A + (Uz[cUID].z - Uz[aUID].z) * B > 0.0)
								{
									A *= -1.0;
									B *= -1.0;
								}
								double d = (A * (Uz[i].r - Uz[aUID].r) + B * (Uz[i].z - Uz[aUID].z)) / sqrt(pow(A, 2.0) + pow(B, 2.0));
								double del = 0.001 * sqrt(pow(A, 2.0) + pow(B, 2.0));
								// if (Uz[i].UID == 60)
								// {
								// 	cout << "Контакт произошел" << "\n";
								// 	cout << d << " " << del << "\n";
								// }
								// if (Uz[i].UID == 4861 && id_ych_gr == 13023)
								// {
								// 	cout << d << " " << del << "\n";
								// }
								if (d >= 0.0 && d <= del)
								{
									if ((fabs(xf1 - Uz[aUID].r) <= 1e-8 && fabs(yf1 - Uz[aUID].z) <= 1e-8) ||
										(fabs(xf1 - Uz[bUID].r) <= 1e-8 && fabs(yf1 - Uz[bUID].z) <= 1e-8))
									{
										if (fabs(xf1 - Uz[aUID].r) <= 1e-8 && fabs(yf1 - Uz[aUID].z) <= 1e-8)
										{
											// if (Uz[i].UID == 60)
											// {
											// 	cout << "Сошелся с узлом a" << "\n";
											// 	cout << xf1 << " " << Uz[aUID].r << " " << xf1 - Uz[aUID].r << "\n";
											// 	cout << yf1 << " " << Uz[aUID].z << " " << yf1 - Uz[aUID].z << "\n";
											// }
											contact(1, Ych, Uz, i, aUID, bUID, lab, del, sinab, cosab, sigmad, id_ych_gr, chislitr, chislitz, omega, gama, beta, Gx, Gy, Rx, Ry, b, flag, Press, time);
										}
										else if (fabs(xf1 - Uz[bUID].r) <= 1e-8 && fabs(yf1 - Uz[bUID].z) <= 1e-8)
										{
											// if (Uz[i].UID == 60)
											// {
											// 	cout << "Сошелся с узлом b" << "\n";
											// 	cout << xf1 << " " << Uz[bUID].r << " " << xf1 - Uz[bUID].r << "\n";
											// 	cout << yf1 << " " << Uz[bUID].z << " " << yf1 - Uz[bUID].z << "\n";
											// }
											contact(-1, Ych, Uz, i, aUID, bUID, lab, del, sinab, cosab, sigmad, id_ych_gr, chislitr, chislitz, omega, gama, beta, Gx, Gy, Rx, Ry, b, flag, Press, time);
										}
									}
									else if (fabs(lab - laf1 - lbf1) <= 1e-8)
									{
										// if (Uz[i].UID == 60)
										// {
										// 	cout << "Зашел посередине" << "\n";
										// 	cout << id_ych_gr + 1 << "\n";
										// }
										contact(1, Ych, Uz, i, aUID, bUID, lab, del, sinab, cosab, sigmad, id_ych_gr, chislitr, chislitz, omega, gama, beta, Gx, Gy, Rx, Ry, b, flag, Press, time);
										// if (Uz[i].UID == 6743 && id_ych_gr == 11657 - 1)
										// {
										// 	cout << G << " " << R << "\n";
										// }
										contact(-1, Ych, Uz, i, aUID, bUID, lab, del, sinab, cosab, sigmad, id_ych_gr, chislitr, chislitz, omega, gama, beta, Gx, Gy, Rx, Ry, b, flag, Press, time);
									}
								}
								else if (d < 0.0)
								{
									double alpha = ((Uz[bUID].r - Uz[aUID].r) * (Uz[i].r - Uz[aUID].r) +
										(Uz[bUID].z - Uz[aUID].z) * (Uz[i].z - Uz[aUID].z)) /
										(pow(Uz[bUID].r - Uz[aUID].r, 2.0) + pow(Uz[bUID].z - Uz[aUID].z, 2.0));
									if (alpha > 0.0 && alpha < 1.0)
									{
										double A = Uz[bUID].z - Uz[aUID].z;
										double B = -(Uz[bUID].r - Uz[aUID].r);
										if ((Uz[cUID].r - Uz[aUID].r) * A + (Uz[cUID].z - Uz[aUID].z) * B > 0.0)
										{
											A *= -1.0;
											B *= -1.0;
										}
										l = A / sqrt(pow(A, 2.0) + pow(B, 2.0));
										m = B / sqrt(pow(A, 2.0) + pow(B, 2.0));
										double Na = l * Uz[aUID].vr + m * Uz[aUID].vz;
										double Nb = l * Uz[bUID].vr + m * Uz[bUID].vz;
										double Nf = l * Uz[i].vr + m * Uz[i].vz;
										double Mb = 0.0;
										double Ma = 0.0;
										double Mf = 0.0;
										for (unsigned int l = 0; l < Uz[bUID].kolvoY.size(); l++)
										{
											Mb += Ych[Uz[bUID].kolvoY[l] - 1].M;
										}
										for (unsigned int l = 0; l < Uz[aUID].kolvoY.size(); l++)
										{
											Ma += Ych[Uz[aUID].kolvoY[l] - 1].M;
										}
										for (unsigned int l = 0; l < Uz[i].kolvoY.size(); l++)
										{
											Mf += Ych[Uz[i].kolvoY[l] - 1].M;
										}
										Mb /= Uz[bUID].kolvoY.size();
										Ma /= Uz[aUID].kolvoY.size();
										Mf /= Uz[i].kolvoY.size();
										double Pl = Ma * Na + Mb * Nb + Mf * Nf;
										double Po = Mb * Nb + alpha * Mf * Nf;
										double Nfp = ((1.0 - alpha) * Mb * Pl - ((1.0 - alpha) * Mb - alpha * Ma) * Po) /
											(Mf * (pow(alpha, 2.0) * Ma + Mb * pow(1.0 - alpha, 2.0)) + Ma * Mb);
										deltaN = Nfp - Nf;
									}
									else
									{
										deltaN = 0;
									}
									flad += 1;
								}
							}
							break;
						}
						case 1:
						case 2:
						{
							aUID = Ych[id_ych_gr].N1ID - 1;
							bUID = Ych[id_ych_gr].N2ID - 1;
							cUID = Ych[id_ych_gr].N3ID - 1;
							double xf1 = 0.0, yf1 = 0.0, lff1 = 0.0;
							normal(Uz, i, aUID, bUID, xf1, yf1, lff1);
							double lab = sqrt(pow(Uz[aUID].r - Uz[bUID].r, 2.0) + pow(Uz[aUID].z - Uz[bUID].z, 2.0));
							double laf1 = sqrt(pow(Uz[aUID].r - xf1, 2.0) + pow(Uz[aUID].z - yf1, 2.0));
							double lbf1 = sqrt(pow(Uz[bUID].r - xf1, 2.0) + pow(Uz[bUID].z - yf1, 2.0));
							// if (Uz[i].UID == 8025 && id_ych_gr == 7789 - 1)
							// {
							// 	cout << "Нашел ячейку с тройным гранич узлом" << "\n";
							// 	cout << aUID << " " << bUID << " " << cUID << "\n";
							// 	cout << lff1 << " " << 1e-4 / lab << "\n";
							// }
							// if (Uz[i].UID == 8025 && id_ych_gr == 7788)
							// {
							// 	cout << aUID + 1 << " " << bUID + 1 << "\n";
							// 	cout << lff1 << " " << 1e-14 / lab << " " << lab - laf1 - lbf1 << "\n";
							// }
							
							if (lff1 <= 1e-8 / lab && fabs(lab - laf1 - lbf1) <= 1e-6)
							{
								// if (Uz[i].UID == 8025)
								// {
								// 	cout << "Контакт с ячейкой: " << id_ych_gr + 1 << "\n";
								// }

								double A = Uz[bUID].z - Uz[aUID].z;
								double B = -(Uz[bUID].r - Uz[aUID].r);
								if ((Uz[cUID].r - Uz[aUID].r) * A + (Uz[cUID].z - Uz[aUID].z) * B > 0.0)
								{
									A *= -1.0;
									B *= -1.0;
								}
								double d = (A * (Uz[i].r - Uz[aUID].r) + B * (Uz[i].z - Uz[aUID].z)) / sqrt(pow(A, 2.0) + pow(B, 2.0));
								double del = 0.001 * sqrt(pow(A, 2.0) + pow(B, 2.0));
								// if (Uz[i].UID == 8025 && id_ych_gr == 7789 - 1)
								// {
								// 	cout << "Проверяем ячейку с тройным гранич узлом" << "\n";
								// 	cout << d << " " << del << "\n";
								// }
								if (d >= 0.0 && d <= del)
								{
									if ((fabs(xf1 - Uz[aUID].r) <= 1e-8 && fabs(yf1 - Uz[aUID].z) <= 1e-8) ||
										(fabs(xf1 - Uz[bUID].r) <= 1e-8 && fabs(yf1 - Uz[bUID].z) <= 1e-8))
									{
										if (fabs(xf1 - Uz[aUID].r) <= 1e-8 && fabs(yf1 - Uz[aUID].z) <= 1e-8)
										{
											contact(1, Ych, Uz, i, aUID, bUID, lab, del, sinab, cosab, sigmad, id_ych_gr, chislitr, chislitz, omega, gama, beta, Gx, Gy, Rx, Ry, b, flag, Press, time);
										}
										else if (fabs(xf1 - Uz[bUID].r) <= 1e-8 && fabs(yf1 - Uz[bUID].z) <= 1e-8)
										{
											contact(-1, Ych, Uz, i, aUID, bUID, lab, del, sinab, cosab, sigmad, id_ych_gr, chislitr, chislitz, omega, gama, beta, Gx, Gy, Rx, Ry, b, flag, Press, time);
										}
									}
									else if (fabs(lab - laf1 - lbf1) <= 1e-8)
									{
										contact(1, Ych, Uz, i, aUID, bUID, lab, del, sinab, cosab, sigmad, id_ych_gr, chislitr, chislitz, omega, gama, beta, Gx, Gy, Rx, Ry, b, flag, Press, time);
										contact(-1, Ych, Uz, i, aUID, bUID, lab, del, sinab, cosab, sigmad, id_ych_gr, chislitr, chislitz, omega, gama, beta, Gx, Gy, Rx, Ry, b, flag, Press, time);
									}
								}
								else if (d < 0.0)
								{
									double alpha = ((Uz[bUID].r - Uz[aUID].r) * (Uz[i].r - Uz[aUID].r) +
										(Uz[bUID].z - Uz[aUID].z) * (Uz[i].z - Uz[aUID].z)) /
										(pow(Uz[bUID].r - Uz[aUID].r, 2.0) + pow(Uz[bUID].z - Uz[aUID].z, 2.0));
									if (alpha > 0.0 && alpha < 1.0)
									{
										double A = Uz[bUID].z - Uz[aUID].z;
										double B = -(Uz[bUID].r - Uz[aUID].r);
										if ((Uz[cUID].r - Uz[aUID].r) * A + (Uz[cUID].z - Uz[aUID].z) * B > 0.0)
										{
											A *= -1.0;
											B *= -1.0;
										}
										l = A / sqrt(pow(A, 2.0) + pow(B, 2.0));
										m = B / sqrt(pow(A, 2.0) + pow(B, 2.0));
										double Na = l * Uz[aUID].vr + m * Uz[aUID].vz;
										double Nb = l * Uz[bUID].vr + m * Uz[bUID].vz;
										double Nf = l * Uz[i].vr + m * Uz[i].vz;
										double Mb = 0.0;
										double Ma = 0.0;
										double Mf = 0.0;
										for (unsigned int l = 0; l < Uz[bUID].kolvoY.size(); l++)
										{
											Mb += Ych[Uz[bUID].kolvoY[l] - 1].M;
										}
										for (unsigned int l = 0; l < Uz[aUID].kolvoY.size(); l++)
										{
											Ma += Ych[Uz[aUID].kolvoY[l] - 1].M;
										}
										for (unsigned int l = 0; l < Uz[i].kolvoY.size(); l++)
										{
											Mf += Ych[Uz[i].kolvoY[l] - 1].M;
										}
										Mb /= Uz[bUID].kolvoY.size();
										Ma /= Uz[aUID].kolvoY.size();
										Mf /= Uz[i].kolvoY.size();
										double Pl = Ma * Na + Mb * Nb + Mf * Nf;
										double Po = Mb * Nb + alpha * Mf * Nf;
										double Nfp = ((1.0 - alpha) * Mb * Pl - ((1.0 - alpha) * Mb - alpha * Ma) * Po) /
											(Mf * (pow(alpha, 2.0) * Ma + Mb * pow(1.0 - alpha, 2.0)) + Ma * Mb);
										deltaN = Nfp - Nf;
									}
									else
									{
										deltaN = 0;
									}
									flad += 1;
								}
							}
							else
							{
								aUID = Ych[id_ych_gr].N2ID - 1;
								bUID = Ych[id_ych_gr].N3ID - 1;
								cUID = Ych[id_ych_gr].N1ID - 1;
								double xf1 = 0.0, yf1 = 0.0, lff1 = 0.0;
								normal(Uz, i, aUID, bUID, xf1, yf1, lff1);
								double lab = sqrt(pow(Uz[aUID].r - Uz[bUID].r, 2.0) + pow(Uz[aUID].z - Uz[bUID].z, 2.0));
								double laf1 = sqrt(pow(Uz[aUID].r - xf1, 2.0) + pow(Uz[aUID].z - yf1, 2.0));
								double lbf1 = sqrt(pow(Uz[bUID].r - xf1, 2.0) + pow(Uz[bUID].z - yf1, 2.0));
								// if (Uz[i].UID == 8025 && id_ych_gr == 7788)
								// {
								// 	cout << aUID + 1 << " " << bUID + 1 << "\n";
								// 	cout << lff1 << " " << 1e-14 / lab << " " << lab - laf1 - lbf1 << "\n";
								// }
								if (lff1 <= 1e-8 / lab && fabs(lab - laf1 - lbf1) <= 1e-6)
								{
									double A = Uz[bUID].z - Uz[aUID].z;
									double B = -(Uz[bUID].r - Uz[aUID].r);
									if ((Uz[cUID].r - Uz[aUID].r) * A + (Uz[cUID].z - Uz[aUID].z) * B > 0.0)
									{
										A *= -1.0;
										B *= -1.0;
									}
									double d = (A * (Uz[i].r - Uz[aUID].r) + B * (Uz[i].z - Uz[aUID].z)) / sqrt(pow(A, 2.0) + pow(B, 2.0));
									double del = 0.001 * sqrt(pow(A, 2.0) + pow(B, 2.0));
									// if (Uz[i].UID == 8025 && id_ych_gr == 7788)
									// {
									// 	cout << d << " " << del << "\n";
									// }
									if (d >= 0.0 && d <= del)
									{
										if ((fabs(xf1 - Uz[aUID].r) <= 1e-8 && fabs(yf1 - Uz[aUID].z) <= 1e-8) ||
											(fabs(xf1 - Uz[bUID].r) <= 1e-8 && fabs(yf1 - Uz[bUID].z) <= 1e-8))
										{
											if (fabs(xf1 - Uz[aUID].r) <= 1e-8 && fabs(yf1 - Uz[aUID].z) <= 1e-8)
											{
												contact(1, Ych, Uz, i, aUID, bUID, lab, del, sinab, cosab, sigmad, id_ych_gr, chislitr, chislitz, omega, gama, beta, Gx, Gy, Rx, Ry, b, flag, Press, time);
											}
											else if (fabs(xf1 - Uz[bUID].r) <= 1e-8 && fabs(yf1 - Uz[bUID].z) <= 1e-8)
											{
												contact(-1, Ych, Uz, i, aUID, bUID, lab, del, sinab, cosab, sigmad, id_ych_gr, chislitr, chislitz, omega, gama, beta, Gx, Gy, Rx, Ry, b, flag, Press, time);
											}
										}
										else if (fabs(lab - laf1 - lbf1) <= 1e-8)
										{
											contact(1, Ych, Uz, i, aUID, bUID, lab, del, sinab, cosab, sigmad, id_ych_gr, chislitr, chislitz, omega, gama, beta, Gx, Gy, Rx, Ry, b, flag, Press, time);
											contact(-1, Ych, Uz, i, aUID, bUID, lab, del, sinab, cosab, sigmad, id_ych_gr, chislitr, chislitz, omega, gama, beta, Gx, Gy, Rx, Ry, b, flag, Press, time);
										}
									}
									else if (d < 0.0)
									{
										double alpha = ((Uz[bUID].r - Uz[aUID].r) * (Uz[i].r - Uz[aUID].r) +
											(Uz[bUID].z - Uz[aUID].z) * (Uz[i].z - Uz[aUID].z)) /
											(pow(Uz[bUID].r - Uz[aUID].r, 2.0) + pow(Uz[bUID].z - Uz[aUID].z, 2.0));
										if (alpha > 0.0 && alpha < 1.0)
										{
											double A = Uz[bUID].z - Uz[aUID].z;
											double B = -(Uz[bUID].r - Uz[aUID].r);
											if ((Uz[cUID].r - Uz[aUID].r) * A + (Uz[cUID].z - Uz[aUID].z) * B > 0.0)
											{
												A *= -1.0;
												B *= -1.0;
											}
											l = A / sqrt(pow(A, 2.0) + pow(B, 2.0));
											m = B / sqrt(pow(A, 2.0) + pow(B, 2.0));
											double Na = l * Uz[aUID].vr + m * Uz[aUID].vz;
											double Nb = l * Uz[bUID].vr + m * Uz[bUID].vz;
											double Nf = l * Uz[i].vr + m * Uz[i].vz;
											double Mb = 0.0;
											double Ma = 0.0;
											double Mf = 0.0;
											for (unsigned int l = 0; l < Uz[bUID].kolvoY.size(); l++)
											{
												Mb += Ych[Uz[bUID].kolvoY[l] - 1].M;
											}
											for (unsigned int l = 0; l < Uz[aUID].kolvoY.size(); l++)
											{
												Ma += Ych[Uz[aUID].kolvoY[l] - 1].M;
											}
											for (unsigned int l = 0; l < Uz[i].kolvoY.size(); l++)
											{
												Mf += Ych[Uz[i].kolvoY[l] - 1].M;
											}
											Mb /= Uz[bUID].kolvoY.size();
											Ma /= Uz[aUID].kolvoY.size();
											Mf /= Uz[i].kolvoY.size();
											double Pl = Ma * Na + Mb * Nb + Mf * Nf;
											double Po = Mb * Nb + alpha * Mf * Nf;
											double Nfp = ((1.0 - alpha) * Mb * Pl - ((1.0 - alpha) * Mb - alpha * Ma) * Po) /
												(Mf * (pow(alpha, 2.0) * Ma + Mb * pow(1.0 - alpha, 2.0)) + Ma * Mb);
											deltaN = Nfp - Nf;
										}
										else
										{
											deltaN = 0;
										}
										flad += 1;
									}

								}
							}
							break;
						}
						default:
							break;
						}
					}

				}

			}
			// if (Uz[i].UID == 3249 || Uz[i].UID == 3019 || Uz[i].UID == 3021)
			// {
			// 	cout << "Узел: " <<  Uz[i].UID << "\n";
			// 	cout << G << " " << R << " " << deltaN << " " << l << " " << m << "\n";
			// }

			if (flad == 0 && flag == 0)
			{
				vrn = Uz[i].vr + deltaT * chislitr / (2.0 * omega) + gama * deltaT;
				vzn = Uz[i].vz + deltaT * chislitz / (2.0 * omega) + beta * deltaT;
				rn = Uz[i].r + vrn * deltaT12;
				zn = Uz[i].z + vzn * deltaT12;
				if (Uz[i].r <= 1e-8)
				{
					vrn = 0.0;
					rn = 0.0;
				}
				
				// if (Uz[i].UID == 7935)
				// {
				// 	cout << "-------\n";
				// 	cout << vrn << " " << Uz[i].vr << "\n";
				// 	cout << vzn << " " << Uz[i].vz <<  "\n";
				// 	cout << rn << " " << Uz[i].r << "\n";
				// 	cout << zn << " " << Uz[i].z << "\n";
				// }
			}
			else if (flag != 0 && flad == 0)
			{
				z = 1 + b / (znamen_m);
				vrn = Uz[i].vr + deltaT * (Gx + Rx);
				vzn = Uz[i].vz + deltaT * (Gy + Ry);
				rn = Uz[i].r + vrn * deltaT12;
				zn = Uz[i].z + vzn * deltaT12;
				//cout << G << " " << R << " " << vrn << " " << vzn << " " << rn << " " << zn << "\n";
				// if (isnan(G) == true || isnan(R) == true)
				// {
				// 	cout << Uz[i].UID << "\n";
				// }
				// if (Uz[i].UID == 7935)
				// {
				// 	cout << "...........\n";
				// 	cout << vrn << " " << Uz[i].vr << "\n";
				// 	cout << vzn << " " << Uz[i].vz <<  "\n";
				// 	cout << rn << " " << Uz[i].r << "\n";
				// 	cout << zn << " " << Uz[i].z << "\n";
				// }
			}
			else if (flad != 0 && flag == 0)
			{
				vrn = Uz[i].vr + deltaT * chislitr / (2.0 * omega) + gama * deltaT + deltaN * l;
				vzn = Uz[i].vz + deltaT * chislitz / (2.0 * omega) + beta * deltaT + deltaN * m;
				rn = Uz[i].r + vrn * deltaT12;
				zn = Uz[i].z + vzn * deltaT12;
				//cout << flad << " " << vrn << " " << vzn << " " << rn << " " << zn << "\n";
				// if (Uz[i].UID == 7935)
				// {
				// 	cout << ",,,,,,,,,\n";
				// 	cout << vrn << " " << Uz[i].vr << "\n";
				// 	cout << vzn << " " << Uz[i].vz <<  "\n";
				// 	cout << rn << " " << Uz[i].r << "\n";
				// 	cout << zn << " " << Uz[i].z << "\n";
				// }
			}
			else if (flag != 0 && flad != 0)
			{
				z = 1 + b / (znamen_m);
				vrn = Uz[i].vr + deltaT * (Gx + Rx) + deltaN * l;
				vzn = Uz[i].vz + deltaT * (Gy + Ry) + deltaN * m;
				rn = Uz[i].r + vrn * deltaT12;
				zn = Uz[i].z + vzn * deltaT12;
			}
			break;
		}
		default:
			break;
		}
	}
}









//===================================
//Скорости деформаций и деформации 
//===================================
void Defor(const Cell& cell, Cell& pseucell, std:: vector<Node>& nodes, const double time) 
{
    //УСМ
    pseucell.A12 = (pseucell.A + cell.A) / 2;
    pseucell.V12 = (pseucell.V + cell.V) / 2;
    if (time != 0.0) 
	{
		pseucell.V12tV = (pseucell.V - cell.V) / (deltaT * pseucell.V12);
    } 
	else 
	{
        pseucell.V12tV = 0;
    }
    
	for (size_t i = 0; i + 1 < cell.NID.size(); ++i)
	{
		Node& node1 = nodes[cell.NID[i] - 1];
		for (size_t j = i + 1; i < cell.NID.size(); ++j)
		{
			Node& node2 = nodes[cell.NID[j] - 1];
			if(i == 0 && j == 2)
			{
				pseucell.err += ApproxNode(node2.vrn, node1.vrn, node2.zn + node2.z, node1.zn + node1.z);
				pseucell.ezz += ApproxNode(node2.vzn, node1.vzn, node1.rn + node1.r, node2.rn + node2.r);
				pseucell.erz += ApproxNode(node2.vzn, node1.vzn, node2.zn + node2.z, node1.zn + node1.z) + 
									ApproxNode(node2.vrn, node1.vrn, node1.rn + node1.r, node2.rn + node2.r);
			}
			else
			{
				pseucell.err += ApproxNode(node1.vrn, node2.vrn, node1.zn + node1.z, node2.zn + node2.z);
				pseucell.ezz += ApproxNode(node2.vzn, node1.vzn, node2.rn + node2.r, node1.rn + node1.r);
				pseucell.erz += ApproxNode(node2.vzn, node1.vzn, node1.zn + node1.z, node2.zn + node2.z) + 
									ApproxNode(node2.vrn, node1.vrn, node2.rn + node2.r, node1.rn + node1.r);
			}
		}
	}
	pseucell.err /= 4 * pseucell.A12;
    pseucell.ezz /= 4 * pseucell.A12;
    pseucell.erz /= 4 * pseucell.A12;
    pseucell.e00 = pseucell.V12tV - pseucell.err - pseucell.ezz;
    
	pseucell.delta_err = pseucell.err * deltaT;
    pseucell.delta_ezz = pseucell.ezz * deltaT;
    pseucell.delta_e00 = pseucell.e00 * deltaT;
    pseucell.delta_erz = pseucell.erz * deltaT;
};

//===================================
//Расчет псевдовязкости
//===================================
void PseudoVis(const Cell& cell, Cell& pseucell) 
{
    if (pseucell.V12tV < 0) 
	{
        double qk = lamdak * lamdak * cell.ro0 * pseucell.A12 * pseucell.V12tV * pseucell.V12tV / pseucell.V12;
        double ql = pseucell.c1 * lamdal * cell.ro0 * std:: sqrt(pseucell.A12) * std:: abs(pseucell.V12tV) / pseucell.V12;
        pseucell.qkomb = qk + ql;
        //Вязкая добавка
        double mu = 0.01 * cell.ro0 * sqrt(pseucell.A12) / pseucell.V12;
        pseucell.qrr = 2 * mu * (pseucell.err - pseucell.V12tV / 3);
        pseucell.qzz = 2 * mu * (pseucell.ezz - pseucell.V12tV / 3);
        pseucell.q00 = 2 * mu * (pseucell.e00 - pseucell.V12tV / 3);
        pseucell.qrz = mu * pseucell.erz;
    } 
	else 
	{
        pseucell.qkomb = 0.0;
        pseucell.qrr = 0.0;
        pseucell.qzz = 0.0;
        pseucell.q00 = 0.0;
        pseucell.qrz = 0.0;
    }
};

//===================================
//Расчет девиаторов напряжения
//===================================
void Stress(const Cell& cell, Cell& pseucell, std:: vector<Node>& nodes, std:: vector<Part>& parts, std:: vector<Mat>& mats)
{
	double deltarz = 0.0;
	for (size_t i = 0; i + 1 < cell.NID.size(); ++i)
	{
		Node& node1 = nodes[cell.NID[i] - 1];
		for (size_t j = i + 1; i < cell.NID.size(); ++j)
		{
			Node& node2 = nodes[cell.NID[j] - 1];
			if(i == 0 && j == 2)
			{
				deltarz += ApproxNode(node2.vzn + node2.vz, node1.vzn + node1.vz, node1.z, node2.z) + 
							ApproxNode(node2.vrn + node2.vr, node1.vrn + node1.vr, node1.r, node2.r);
			}
			else
			{
				deltarz += ApproxNode(node2.vzn + node2.vz, node1.vzn + node1.vz, node2.z, node1.z) + 
							ApproxNode(node2.vrn + node2.vr, node1.vrn + node1.vr, node2.r, node1.r);
			}
		}
	}
    deltarz /= 4 * cell.A;
    
	double drr = cell.Drr * deltarz * deltaT; 
    double drz = 0.5 * (cell.Dzz - cell.Drr) * deltarz * deltaT;
	double& G = mats[parts[cell.PID - 1].n_mat - 1].G;
	double& sigmay = mats[parts[cell.PID - 1].n_mat - 1].sigmay;
    
	double intermval = - pseucell.V12tV * deltaT / 3.0;
	pseucell.Drr = cell.Drr + 2.0 * G * (pseucell.delta_err + intermval) + drr;
    pseucell.Dzz = cell.Dzz + 2.0 * G * (pseucell.delta_ezz + intermval) - drr;
    pseucell.D00 = cell.D00 + 2.0 * G * (pseucell.delta_e00 + intermval);
    pseucell.Drz = cell.Drz + G * pseucell.delta_erz + drz;
    
    double f = 2 * (pseucell.Drr * pseucell.Drr + pseucell.Dzz * pseucell.Dzz + pseucell.Drz * pseucell.Drz + pseucell.Drr * pseucell.Dzz);
    if (f > 2 * sigmay * sigmay / 3) {
        double F = sigmay * sqrt(2 / (3 * f));
        pseucell.Drr *= F;
        pseucell.Dzz *= F;
        pseucell.Drz *= F;
        pseucell.D00 *= F;
    }

    //Напряжения 
    pseucell.sigma_rr = (pseucell.Drr + pseucell.qrr) - (pseucell.p + pseucell.qkomb);
    pseucell.sigma_zz = (pseucell.Dzz + pseucell.qzz) - (pseucell.p + pseucell.qkomb);
    pseucell.sigma_00 = (pseucell.D00 + pseucell.q00) - (pseucell.p + pseucell.qkomb);
    pseucell.sigma_rz = (pseucell.Drz + pseucell.qrz); 
	
};

//===================================
//Уравнение энергии
//===================================
void Energy(Cell& cell, Cell& pseucell) 
{
    double D12rr = 0.5 * (pseucell.Drr + cell.Drr); 
    double D12zz = 0.5 * (pseucell.Dzz + cell.Dzz); 
    double D1200 = 0.5 * (pseucell.D00 + cell.D00); 
    double D12rz = 0.5 * (pseucell.Drz + cell.Drz); 
    double Z = pseucell.V12 * (D12rr * pseucell.err + D12zz * pseucell.ezz + D1200 * pseucell.e00 + D12rz * pseucell.erz);
    pseucell.E = cell.E - (cell.p / 2.0 + pseucell.qkomb) * (pseucell.V - cell.V) + Z * deltaT;
}

