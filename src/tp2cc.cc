// On resout l'equation de la chaleur (de Laplace) sur un domaine regulier.
#include <iostream>
#include <iterator>
#include <fstream>
#include "Array2D.hpp"
#include <mpi.h>
#include <math.h>
#include<vector>
#include <string>

// Sauvegarde d'une matrice dans un fichier texte
void save(Array2D<double> &matrix, std::string name) {
	std::ofstream file(name.c_str());
	for (int iY=0; iY<matrix.sizeY(); ++iY) {
		copy(&matrix.data()[iY*matrix.sizeX()], &matrix.data()[iY*matrix.sizeX()]+matrix.sizeX(),
		std::ostream_iterator<double>(file, " "));
		file << "\n";
	}
}

void display(std::vector<double> a){
	std::cout << "printing vector" << std::endl;
		for(int i = 0; i < a.size(); i++){
			std::cout << a.at(i) << "-";
		}
		std::cout << std::endl;
}

void displayInt(std::vector<int> a){
	std::cout << "printing vector" << std::endl;
		for(int i = 0; i < a.size(); i++){
			std::cout << a.at(i) << "-";
		}
		std::cout << std::endl;
}
void displayWithPointer(double *a, int x){
	std::cout << "printing vector with pointer" << std::endl;
		for(int i = 0; i < x; i++){
			std::cout << *a++ << ",";
		}
		std::cout << std::endl;
}
void displayMat(Array2D<double> mat, int dimX, int dimY){
	std::cout << "printing whole matrix" << std::endl;
	for(int iY = 0; iY < dimY; iY++){
		for(int iX = 0;iX < dimX; iX++){
			std::cout << mat(iX,iY) << ",";
		}
		std::cout << std::endl;
	}
		std::cout << std::endl;
}

void solveNormally(Array2D<double> heat,Array2D<double> tmp, int maxT, int dimX, int dimY, std::string fileName){

	  for (int iX=0; iX<dimX; iX++) {      // conditions aux bords:
	      heat(iX,0) = 0;                 // 0 en haut
	      heat(iX,dimY-1) = 1;            // 1 en bas
	      tmp(iX,0) = 0;                  // 0 en haut
	      tmp(iX,dimY-1) = 1;             // 1 en bas
	  }
	  for (int iY=0; iY<dimY; iY++) {
	      heat(0,iY)      = 0.;           // 0 a gauche
	      heat(dimX-1,iY) = 1.;           // 1 a droite
	      tmp(0,iY)      = 0.;            // 0 a gauche
	      tmp(dimX-1,iY) = 1.;            // 1 a droite
	  }

	  for (int iT=0; iT<maxT; iT++) {      // boucle principale : on fait maxT iterations
	    for (int iY=1; iY<dimY-1; iY++) {  // on itere a l'interieur du domaine
	      for (int iX=1; iX<dimX-1; iX++) {
		//std::cout<< "up: " << heat(iX,iY-1) << "down: " << heat(iX,iY+1) << "left: " << heat(iX-1,iY) << "right: " <<heat(iX+1,iY) <<std::endl;
	        tmp(iX,iY) = 0.25*( heat(iX-1,iY) + heat(iX+1,iY)+
	                            heat(iX,iY-1) + heat(iX,iY+1) );
	      }
	    }
	    heat.unsafeSwap(tmp);              // Les deux matrices sont interverties (ceci ne fait que copier
	                                       // deux pointeurs, ca ne coute donc pas trop de temps)
	  }

	  save(heat, fileName);
}
int main(int argc, char **argv) {


	int myRank, nproc;
	MPI_Init(&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);

	//std::cout << "argc: "<< argc << std::endl;

	const int dimX = atoi(argv[1]);
	const int dimY = atoi(argv[2]);
	const int maxT = atoi(argv[3]);
	//const char* ind = argv[4];
	//std::cout << "maxT: " << maxT << std::endl;

	if(nproc > dimY){
		nproc = dimY;
		//return 0;
		MPI_Finalize();
		std::cout << "number of processors is larger than the number of lines which is not valid: " << std::endl;
		return 0;
	}
	int domainSize = ceil((double)dimY/nproc);
	//std::cout << "dimY: "<<dimY<< "domainSize: "<<domainSize<< ceil(8/3) <<std::endl;
	int lastRank = ceil((double)dimY/domainSize);
	lastRank = lastRank - 1;


	Array2D<double> heat(dimX, domainSize,0); // La matrice de la chaleur
	Array2D<double> tmp(dimX, domainSize,-1);  // Une matrice temporaire
	Array2D<double> finalData(1, 1,-1);



	std::cout <<"rank: "<< myRank << "dimY: " << dimY << " dimX: "<<dimX<< " domainSize: "<<domainSize<< std::endl;


	std::vector<int> displacements;
	std::vector<int> recvCounts;
	int disp = 0;
	//initializing displacements and recvCounts
	for(int i=0; i <= lastRank; i++){

		if(i < lastRank){
			recvCounts.push_back(domainSize*dimX);
			if(i == 0){
				displacements.push_back(0);
			}
			else{
				disp = disp + recvCounts[i - 1];
				displacements.push_back(disp);
			}
		}
		else if(i == lastRank){
			int domain = dimY - domainSize*myRank;
			recvCounts.push_back(domain*dimX);
			disp = disp + recvCounts[i - 1];
			displacements.push_back(disp);
		}
	}
	//std::cout <<"rank: "<< myRank << "do you get here2??" << std::endl;

	//if rank is 0

	std::string fileName = "chaleur.dat";
	//fileName = fileName + "_";
	//fileName.append(ind);
	//fileName = fileName + ".dat";
	if(myRank == 0){

		if(nproc == 1){

			//std::string fileNameOrig = "chaleur.dat";
			//fileNameOrig = fileNameOrig + "_";
			//fileNameOrig.append(ind);
			//fileNameOrig = fileNameOrig + ".dat";

			solveNormally(heat, tmp, maxT, dimX, dimY, fileName);
		}
		else{
			finalData.resize(dimX, dimY);
			//std::cout <<"rank: "<< myRank << "I'm the zero rank" << std::endl;
			//std::cout <<"rank: "<< myRank << "do you get here2.1??" << std::endl;
			double right, left, up, down;
			for (int iX=0; iX<dimX; iX++) {      // conditions aux bords:
				heat(iX,0) = 0;                 // 0 en haut
				tmp(iX,0) = 0;                  // 0 en haut
			}
			//std::cout <<"rank: "<< myRank << "do you get here2.2??" << std::endl;
			for (int iY=0; iY<domainSize; iY++){
				//std::cout <<"rank: "<< myRank <<"dimX: " << dimX << "iY: " << iY << "do you get here2.22??" << std::endl;
				heat(0,iY) = 0.;           // 0 a gauche
				tmp(0,iY) = 0.;            // 0 a gauche
			}
			//std::cout <<"rank: "<< myRank << "do you get here2.3??" << std::endl;
			for (int iY=0; iY<domainSize; iY++){
				//std::cout <<"rank: "<< myRank <<"dimX: " << dimX << "iY: " << iY << "do you get here2.4??" << std::endl;
				heat(dimX-1,iY) = 1.;           // 1 a droite
				tmp(dimX-1,iY) = 1.;            // 1 a droite
			}

			//save(heat, "initheat.dat");


			//std::cout <<"rank: "<< myRank << "do you get here2.5??" << std::endl;
			std::vector<double> below(dimX,-1);
			for (int iT=0; iT<maxT; iT++) {      // boucle principale : on fait maxT iterations
				//recv() la ligne en dessous
				//send() la ligne  tout en bas
				double* sendBuffer = heat.data()+ (domainSize - 1)*dimX;// I wonder it the format is correct
				//std::vector<double> sendBuffer(100);
				int sendCount = dimX;
				int destRank = 1;

				//std::cout <<"rank: "<< myRank << "do you get here3??" << std::endl;

				int recvCount = dimX;
				int source = 0;
				MPI_Status status1, status2;
				MPI_Request request1, request2;

				//std::cout <<"rank: "<< myRank << "do you get here4??" << std::endl;

				//std::cout <<"iter: "<< iT <<"rank: "<< myRank << "pointer on sendbuffer: "<< sendBuffer<< std::endl;
				//std::cout <<"iter: "<< iT <<"rank: "<< myRank << "pointer on above: "<< below.data() << std::endl;
				//std::cout <<"iter: "<< iT <<"rank: "<< myRank << "did you send below?" << "destRank: "<< destRank << std::endl;
				//MPI_Sendrecv(sendBuffer, sendCount, MPI_DOUBLE, destRank, 0, below.data(), recvCount, MPI_DOUBLE, source, 0, MPI_COMM_WORLD, status);
				//std::cout << "reading the values of sendBuffer of rank 0" << std::endl;
				//displayWithPointer(sendBuffer, dimX);
				//std::cout << "reading the values of heat.data() of rank 0" << std::endl;
				//displayWithPointer(heat.data(), dimX);

					std::cout <<"************************* iter: "<<iT<<" rank: "<< myRank << "before send and receive" << "destRank: "<< destRank << std::endl;
					int sent = MPI_Isend(sendBuffer, sendCount, MPI_DOUBLE, destRank, 0, MPI_COMM_WORLD, &request1);
					int received = MPI_Irecv( below.data(), recvCount, MPI_DOUBLE, destRank, 0, MPI_COMM_WORLD, &request2);
					//MPI_Barrier(MPI_COMM_WORLD);
					MPI_Wait(&request1, &status1);
					MPI_Wait(&request2, &status2);
					std::cout <<"************************* iter: "<<iT<<" rank: "<< myRank << "after send and receive" << "destRank: "<< destRank << std::endl;
					//std::cout <<"iter: "<< iT <<"rank: "<< myRank << "sent in what situation?" << sent << std::endl;
					//std::cout <<"iter: "<< iT <<"rank: "<< myRank << "received in what situation?" << received<< std::endl;
					//std::cout <<"iter: "<< iT <<"rank: "<< myRank << "sent the below" << std::endl;
					//std::cout <<"iter: "<< iT <<"rank: "<< myRank << "below value? " << below.at(0) << "doma: " << domainSize << std::endl;

				for (int iY= 1; iY < domainSize; iY++) {  // on itere a l'interieur du domaine
					for (int iX=1; iX<dimX-1; iX++) {
						//std::cout <<"iter: "<< iT <<"rank: "<< myRank << "problem here? " << std::endl;
						left = heat(iX-1,iY);
						right = heat(iX+1,iY);
						up = heat(iX,iY-1);

						if(iY == domainSize-1){
							down = below.at(iX);
						}
						else{
							down = heat(iX,iY+1);
						}
						//std::cout<< "up: " << up << "down: " << down << "left: " << left << "right: " << right << std::endl;
						tmp(iX,iY) = 0.25*( left + right + down + up);
						//std::cout << "iter: "<< iT <<" rank: "<< myRank << "result" << tmp(iX,iY) << std::endl;
					}
				}

				//std::cout <<"rank: "<< myRank << "did one of the iterations" << std::endl;
				heat.unsafeSwap(tmp);
				//displayMat(heat, dimX, domainSize);               // Les deux matrices sont interverties (ceci ne fait que copier
				// deux pointeurs, ca ne coute donc pas trop de temps)

			}
			std::cout <<"*******************************rank: "<< myRank << "before gather" << std::endl;
			//MPI_Barrier(MPI_COMM_WORLD);

			MPI_Gatherv(heat.data(), domainSize*dimX, MPI_DOUBLE, finalData.data(), recvCounts.data(), displacements.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
			//MPI_Barrier(MPI_COMM_WORLD);
			//std::cout << "recvCounts: "<< std::endl;
			//displayInt(recvCounts);

			//std::cout << "displacements: "<< std::endl;
			//displayInt(displacements);

			std::cout <<"*******************************rank: "<< myRank << "after gather" << std::endl;
			save(finalData, fileName);

		}


	}
	// if rank in between
	else if(myRank < lastRank){
		//std::cout <<"rank: "<< myRank << "I'm in the between ranks" << std::endl;
		double right, left, up, down;
		for (int iY=0; iY<domainSize; iY++){
			heat(0,iY)      = 0.;           // 0 a gauche
			tmp(0,iY)      = 0.;            // 0 a gauche
		}
		for (int iY=0; iY<domainSize; iY++){
			heat(dimX-1,iY) = 1.;           // 1 a droite
			tmp(dimX-1,iY) = 1.;            // 1 a droite
		}

		//std::cout <<"rank: "<< myRank << "do you get before iterations??" << std::endl;
		std::vector<double> above(dimX,-1);
		std::vector<double> below(dimX,-1);
		for (int iT=0; iT<maxT; iT++) {      // boucle principale : on fait maxT iterations

			//send sa ligne en haut
			//recv ligne au dessus
			double* sendBuffer = heat.data();// I wonder it the format is correct
			//std::vector<double> sendBuffer(100);
			int sendCount = dimX;
			int destRank = myRank - 1;


			int recvCount = dimX;
			int source = myRank;
			MPI_Status  status1, status2, status3, status4;
			MPI_Request request1, request2;


			//std::cout <<"iter: "<< iT <<"rank: "<< myRank << "pointer on sendbuffer: "<< sendBuffer<< std::endl;
			//std::cout <<"iter: "<< iT <<"rank: "<< myRank << "pointer on above: "<< above.data() << std::endl;
			//std::cout <<"iter: "<< iT <<"rank: "<< myRank << "did you send above?" << "destRank: "<< destRank << std::endl;
			//MPI_Sendrecv(sendBuffer, sendCount, MPI_DOUBLE, destRank, 0, above.data(), recvCount, MPI_DOUBLE, source, 0, MPI_COMM_WORLD, status);
			//std::cout << "before reading the values of above of rank 1" << std::endl;
			//displayWithPointer(above.data(), dimX);
			std::cout <<"************************* iter: "<<iT<<" rank: "<< myRank << "before send and receive" << "destRank: "<< destRank << std::endl;
			int sent = MPI_Isend(heat.data(), sendCount, MPI_DOUBLE, destRank, 0, MPI_COMM_WORLD, &request1);
			int received = MPI_Irecv( above.data(), recvCount, MPI_DOUBLE, destRank, 0, MPI_COMM_WORLD, &request2);
			MPI_Wait(&request1, &status1);
			MPI_Wait(&request2, &status2);
			//MPI_Barrier(MPI_COMM_WORLD);
			std::cout <<"************************* iter: "<<iT<<" rank: "<< myRank << "after send and receive" << "destRank: "<< destRank << std::endl;
			//std::cout << "after reading the values of above of rank 1" << std::endl;
			//displayWith(above.data(), dimX);
			//std::cout << "after reading the values of above of rank 1 without pointer" << std::endl;
			//display(above);
			// std::cout <<"iter: "<< iT <<"rank: "<< myRank << "sent in what situation?" << sent << std::endl;
			// std::cout <<"iter: "<< iT <<"rank: "<< myRank << "received in what situation?" << received<< std::endl;
			// std::cout <<"iter: "<< iT <<"rank: "<< myRank << "above? " << above.at(2)<< std::endl;
			// std::cout <<"iter: "<< iT <<"rank: "<< myRank << "sent above" << std::endl;

			//recv ligne en dessous
			//send sa ligne en bas

			sendCount = dimX;
			destRank = myRank + 1;

			recvCount = dimX;
			source = myRank;

			std::cout <<"************************* iter: "<<iT<<" rank: "<< myRank << "before send and receive" << "destRank: "<< destRank << std::endl;
			//MPI_Sendrecv(sendBuffer, sendCount, MPI_DOUBLE, destRank, 0, below.data(), recvCount, MPI_DOUBLE, source, 0, MPI_COMM_WORLD, status);
			sent = MPI_Isend(heat.data()+ (domainSize-1)*dimX, sendCount, MPI_DOUBLE, destRank, 0, MPI_COMM_WORLD, &request1);
			received = MPI_Irecv( below.data(), recvCount, MPI_DOUBLE, destRank, 0, MPI_COMM_WORLD, &request2);
			MPI_Wait(&request1, &status3);
			MPI_Wait(&request2, &status4);
			std::cout <<"************************* iter: "<<iT<<" rank: "<< myRank << "after send and receive" << "destRank: "<< destRank << std::endl;
			//display(below);
			//std::cout <<"rank: "<< myRank << "sent below" << std::endl;
			for (int iY= 0; iY < domainSize; iY++) {  // on itere a l'interieur du domaine
				for (int iX=1; iX<dimX-1; iX++) {
					//std::cout <<"rank: "<< myRank << "here0" << std::endl;
					left = heat(iX-1,iY);
					right = heat(iX+1,iY);
					//std::cout <<"rank: "<< myRank << "here0.5" << std::endl;
					//std::cout <<"rank: "<< myRank << "here0.5" <<"doma: "<< domainSize <<"iY: "<< iY << std::endl;
					if(iY == 0 && iY == domainSize-1){
						down = below.at(iX);
						//std::cout <<"rank: "<< myRank << "here0.6" << std::endl;
						up = above.at(iX);
						//std::cout <<"rank: "<< myRank << "here0.65" << std::endl;
					}
					else if(iY == domainSize-1){
						down = below.at(iX);
						//std::cout <<"rank: "<< myRank << "here0.7" << std::endl;
						up = heat(iX,iY-1);
						//std::cout <<"rank: "<< myRank << "here1" << std::endl;
					}
					else if(iY == 0){
						down = heat(iX,iY+1);
						up = above.at(iX);
						//std::cout <<"rank: "<< myRank << "here2" << std::endl;
					}
					else{
						up = heat(iX,iY-1);
						down = heat(iX,iY+1);
						//std::cout <<"rank: "<< myRank << "here3" << std::endl;
					}
					//std::cout<< "up: " << up << "down: " << down << "left: " << left << "right: " << right << std::endl;
					tmp(iX,iY) = 0.25*( left + right + down + up);
					//std::cout << "iter: "<< iT << " rank: "<< myRank << "result" << tmp(iX,iY) << std::endl;
				}
			}
			heat.unsafeSwap(tmp);              // Les deux matrices sont interverties (ceci ne fait que copier
			//displayMat(heat, dimX, domainSize);
			// deux pointeurs, ca ne coute donc pas trop de temps)
		}
		std::cout <<"*******************************rank: "<< myRank << "before gather" << std::endl;
		//MPI_Barrier(MPI_COMM_WORLD);
		MPI_Gatherv(heat.data(), domainSize*dimX, MPI_DOUBLE, finalData.data(), recvCounts.data(), displacements.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
		std::cout <<"*******************************rank: "<< myRank << "after gather" << std::endl;


	}
	// lastRank
	else if(myRank == lastRank){
		//std::cout <<"rank: "<< myRank << "I'm in the last rank" << std::endl;
		double right, left, up, down;
		int domain = dimY - domainSize*myRank;
		if(domainSize > domain){
				domainSize = domain;
		}

		for (int iX=0; iX<dimX; iX++) {      // conditions aux bords:
			heat(iX,domainSize-1) = 1;             // 1 en bas
			tmp(iX,domainSize-1) = 1;             // 1 en bas
		}
		for (int iY=0; iY<domainSize; iY++){
			heat(0,iY) = 0.;           // 0 a gauche
			tmp(0,iY) = 0.;            // 0 a gauche
		}
		for (int iY=0; iY<domainSize; iY++){
			heat(dimX-1,iY) = 1.;           // 1 a droite
			tmp(dimX-1,iY) = 1.;            // 1 a droite
		}

		//std::cout <<"rank: "<< myRank << "did you get before iteration?" << std::endl;
		std::vector<double> above(dimX,-1);
		for (int iT=0; iT<maxT; iT++) {      // boucle principale : on fait maxT iterations
			//send la ligne tout en haut
			//recv la ligne au dessus

			int sendCount = dimX;
			int destRank = myRank - 1;

			int recvCount = dimX;
			int source = myRank;
			MPI_Status status1, status2;
			MPI_Request request1, request2;

			//std::cout <<"iter: "<< iT <<"rank: "<< myRank << "pointer on above: "<< above.data() << std::endl;
			//std::cout <<"iter: "<< iT <<"rank: "<< myRank << "dimX: "<< dimX << std::endl;
			std::cout <<"************************* iter: "<< iT <<" rank: "<< myRank << "before send and recv" << "destRank: "<< destRank << std::endl;
			//MPI_Sendrecv(sendBuffer, sendCount, MPI_DOUBLE, destRank, 0, above.data(), recvCount, MPI_DOUBLE, source, 0, MPI_COMM_WORLD, status);
			int sent = MPI_Isend(heat.data(), sendCount, MPI_DOUBLE, destRank, 0, MPI_COMM_WORLD, &request1);
			int received = MPI_Irecv( above.data(), recvCount, MPI_DOUBLE, destRank, 0, MPI_COMM_WORLD, &request2);
			MPI_Wait(&request1, &status1);
			MPI_Wait(&request2, &status2);
			std::cout <<"************************* iter: "<< iT <<" rank: "<< myRank << "after send and recv" << "destRank: "<< destRank  << std::endl;
			//std::cout <<"iter: "<< iT <<"rank: "<< myRank << "received in what situation? " << received<< std::endl;
			//std::cout <<"iter: "<< iT <<"rank: "<< myRank << "above? " << above.at(2)<< std::endl;
			//std::cout <<"iter: "<< iT <<"rank: "<< myRank << "sent above" << std::endl;


			for (int iY= 0; iY < domainSize-1; iY++){  // on itere a l'interieur du domaine
				for (int iX=1; iX<dimX-1; iX++){

					left = heat(iX-1,iY);
					//std::cout <<"rank: "<< myRank << "here1" << std::endl;
					right = heat(iX+1,iY);
					//std::cout <<"rank: "<< myRank << "here2" << std::endl;
					down = heat(iX,iY+1);
					//std::cout <<"rank: "<< myRank << "here3" << std::endl;

					if(iY == 0){
						//std::cout <<"rank: "<< myRank << "here3.2" << std::endl;
						up = above.at(iX);
						//std::cout <<"rank: "<< myRank << "here3.5" << std::endl;
					}
					else{
						//std::cout <<"rank: "<< myRank << "here3.7" << std::endl;
						up = heat(iX,iY-1);
						//std::cout <<"rank: "<< myRank << "here3.8" << std::endl;
					}
					//std::cout <<"rank: "<< myRank << "here4" << std::endl;
					//std::cout<< "up: " << up << "down: " << down << "left: " << left << "right: " << right << std::endl;
					tmp(iX,iY) = 0.25*( left + right + down + up);
					//std::cout << "iter: "<< iT << " rank: "<< myRank << "result" << tmp(iX,iY) << std::endl;
				}
			}
			heat.unsafeSwap(tmp);
			//displayMat(heat, dimX, domainSize);             // Les deux matrices sont interverties (ceci ne fait que copier
			//std::cout <<"rank: "<< myRank << "here4" << std::endl;
			// deux pointeurs, ca ne coute donc pas trop de temps)
		}
		std::cout <<"*******************************rank: "<< myRank << "before gather" << std::endl;
		//MPI_Barrier(MPI_COMM_WORLD);
		MPI_Gatherv(heat.data(), domainSize*dimX, MPI_DOUBLE, finalData.data(), recvCounts.data(), displacements.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
		std::cout <<"*******************************rank: "<< myRank << "after gather" << std::endl;
		//std::cout <<"*******************************rank: "<< myRank << "after gather" << std::endl;


	}

		MPI_Finalize();

	}
