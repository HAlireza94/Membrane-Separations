#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
using namespace std;
template<typename T>


vector<double> linspace(T start_in, T end_in, int num_in)
{

    vector<double> linspaced;

    double start = static_cast<double>(start_in);
    double end = static_cast<double>(end_in);
    double num = static_cast<double>(num_in);

    if (num == 0) { return linspaced; }
    if (num == 1)
    {
        linspaced.push_back(start);
        return linspaced;
    }

    double delta = (end - start) / (num - 1);

    for(int i=0; i < num-1; ++i)
    {
        linspaced.push_back(start + delta * i);
    }
    linspaced.push_back(end); // I want to ensure that start and end
    // are exactly the same as the input
    return linspaced;
}

void print_vector(vector<double> vec)
{
    cout << "size: " << vec.size() << endl;
    for (double d : vec)
        cout << d << " ";
    cout << std::endl;
}

void print2DVector(const vector<vector<double>>& vec) {
    for (const auto& row : vec) {
        for (const auto& elem : row) {
            cout << elem << " ";
        }
        cout << std::endl; // End of row
    }
}

int main()
{
    // Simulation Details
    int const Lx = 20;// nm
    float const Ly = 6.39; // nm
    float const Lz = 4.919; // # nm
    int const numberofWater = 8951;
    int const numberofAtoms = 4 * numberofWater;
    int const numberofSteps = 25000000;
    float const dt = 0.002;
    int const productionTime = numberofSteps * dt * 0.001; // ns
    int const stepXtc = 200; // step size to save on xtc
    int const stepeachframeTime = stepXtc * dt; //ps
    //char const targetAtom = 'OW';
    //int numberofFrames = (numberofSteps / stepXtc) + 1;

    vector<double> binSize = linspace(0, 20, 81);
    float const delta = binSize[1] - binSize[0];



    string line;
    ifstream infile("working.gro");

    if (!infile.is_open()) {
        std::cerr << "Error opening file" << std::endl;
        return 1;
    }

    vector<float> X1;
    while (getline(infile, line)) {
        istringstream iss(line);
        vector<string> tokens;
        string token;

        while (iss >> token) {
            tokens.push_back(token);
        }

        if (tokens.size() == 6) {
            X1.push_back(stof(tokens[3]));
        } else if (tokens.size() == 5) {
            X1.push_back(stof(tokens[2]));
        }
    }

    infile.close();




    int numberofFrames = (X1.size() / numberofAtoms) + 1;

    vector<vector<double>> densityNumber(binSize.size() - 1,
                                         vector<double>(numberofFrames - 1, 0.0));

    for (int i = 1; i < numberofFrames; ++i) {
        for (size_t j = 0; j < binSize.size() - 1; ++j) {
            int p = 0;
            for (int k = i * numberofAtoms - numberofAtoms; k < i * numberofAtoms; k += 4) {
                if (binSize[j] < X1[k] && X1[k] < binSize[j + 1]) {
                    ++p;
                }
            }
            densityNumber[j][i - 1] = static_cast<double>(p) / (delta * Ly * Lz);
        }
    }

    vector<double> rowSums(densityNumber.size(), 0.0);
    for (size_t i = 0; i < densityNumber.size(); ++i) {
        double sum = 0.0;
        for (size_t j = 0; j < densityNumber[i].size(); ++j) {
            sum += densityNumber[i][j];
        }
        rowSums[i] = sum/(numberofFrames-1);
    }

    ofstream output("Density.xvg");
    // Optional: code to print or use the row sums

    for (int i = 0; i < binSize.size() -1; i++){

        output << binSize[i] << "   " << rowSums[i] << '\n';
    }

    output.close();


    




    return 0;
}