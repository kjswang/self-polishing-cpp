#ifndef _LOADPARAMETER_H_
#define _LOADPARAMETER_H_

// basic file operations

/*
* Parameter loading library 
* Modifier - 
* Date - 
*/

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "structure.h"
#include <map>
#include "math_pack.h"
#include "yaml-cpp/yaml.h"  // IWYU pragma: keep
#include <armadillo>
// #include <typeinfo>
// #include <Eigen/Dense>

// using namespace Eigen;
using namespace std;
using namespace arma;


void setMatrix(MatrixXd& row, string& data){
    /*
    * load matrix row from txt file
    * data is a matrix row (string type)
    * input data format: a,b,c,..,z
    * abc...z are row elements
    */
    // words count 
    int cnt = 0;
    vector<float> row_tmp;
    string element = ""; 
    for (auto x : data){
        if (x == ',') 
       {    
            // get row element
            // cout << element << endl;
            row_tmp.push_back(stof(element)); // push back next element in float
            element = ""; // empty element
            cnt++;
       } 
       else
       { 
           element = element + x; 
       } 
    }
    // the last element doesn't have ","" behind it
    row_tmp.push_back(stof(element)); // push back next element in float
    // initilize matrix 
    row.resize(1,row_tmp.size());
    // push row vector value to Eigen::Matrix from input
    cnt = 0;
    for (vector<float>::const_iterator i = row_tmp.begin(); i != row_tmp.end(); ++i){
        row(0,cnt) = double(*i); // row is a one row matrix
        cnt++;
    }
}

void setMatrix_arma(mat& row, string& data){
    int cnt = 0;
    vector<float> row_tmp;
    string element = ""; 
    for (auto x : data){
        if (x == ','){    
            row_tmp.push_back(stof(element));
            element = "";
            cnt++;
        } 
        else {
            element = element + x; 
        } 
    }
    row_tmp.push_back(stof(element));
    row.resize(1, row_tmp.size());
    cnt = 0;
    for (vector<float>::const_iterator i = row_tmp.begin(); i != row_tmp.end(); ++i){
        row(0, cnt) = double(*i);
        cnt++;
    }
}

void loadWeldLeft(MatrixXd& weldleft){
    fstream file("parameter/wp_setup/weld_left.txt");
    string data;
    inputfield line;
    map<string, vector<float> > dicts;
    MatrixXd row_tmp;
    MatrixXd weld_left_tmp;
    int cnt =0;
    while (std::getline(file,data)) {
        // skip the empty line 
        if (data.length()==0)
            continue;
        // customized matrix setting
        // get the matrix elements in string
        setMatrix(row_tmp, data);
        // concantenate current row with old weld traj
        if (cnt == 0){
            weld_left_tmp = row_tmp;
            cnt++;
        } // the first row 
        else {
            weld_left_tmp = Vcat(weld_left_tmp, row_tmp);
            cnt++;
        }
    }
    weldleft = weld_left_tmp;
}

void loadWeldLeft_arma(mat& PC_idx){
    ifstream file("parameter/wp_setup/weld_left.txt");
    string data;
    mat row_tmp;
    mat PC_idx_tmp;
    int cnt = 0;
    while (getline(file, data)) {
        if (data.length() == 0)
            continue;
        setMatrix_arma(row_tmp, data);
        if (cnt == 0){
            PC_idx_tmp = row_tmp;
            cnt++;
        } else {
            PC_idx_tmp = join_cols(PC_idx_tmp, row_tmp);
            cnt++;
        }
    }
    PC_idx = PC_idx_tmp;
}

void loadWeldRight(MatrixXd& weldright){
    fstream file("parameter/wp_setup/weld_right.txt");
    string data;
    inputfield line;
    map<string, vector<float> > dicts;
    MatrixXd row_tmp;
    MatrixXd weld_right_tmp;
    int cnt =0;
    while (std::getline(file,data)) {
        // skip the empty line 
        if (data.length()==0)
            continue;
        // customized matrix setting
        // get the matrix elements in string
        setMatrix(row_tmp, data);
        // concantenate current row with old weld traj
        if (cnt == 0){
            weld_right_tmp = row_tmp;
            cnt++;
        } // the first row 
        else {
            weld_right_tmp = Vcat(weld_right_tmp, row_tmp);
            cnt++;
        }
    }
    weldright = weld_right_tmp;
}

void loadWeldRight_arma(mat& PC_idx){
    ifstream file("parameter/wp_setup/weld_right.txt");
    string data;
    mat row_tmp;
    mat PC_idx_tmp;
    int cnt = 0;
    while (getline(file, data)) {
        if (data.length() == 0)
            continue;
        setMatrix_arma(row_tmp, data);
        if (cnt == 0){
            PC_idx_tmp = row_tmp;
            cnt++;
        } else {
            PC_idx_tmp = join_cols(PC_idx_tmp, row_tmp);
            cnt++;
        }
    }
    PC_idx = PC_idx_tmp;
}

void loadWeldTop(MatrixXd& weldtop){
    fstream file("parameter/wp_setup/weld_top.txt");
    string data;
    inputfield line;
    map<string, vector<float> > dicts;
    MatrixXd row_tmp;
    MatrixXd weld_top_tmp;
    int cnt =0;
    while (std::getline(file,data)) {
        // skip the empty line 
        if (data.length()==0)
            continue;
        // customized matrix setting
        // get the matrix elements in string
        setMatrix(row_tmp, data);
        // concantenate current row with old weld traj
        if (cnt == 0){
            weld_top_tmp = row_tmp;
            cnt++;
        } // the first row 
        else {
            weld_top_tmp = Vcat(weld_top_tmp, row_tmp);
            cnt++;
        }
    }
    weldtop = weld_top_tmp;
}

void loadWeldTop_arma(mat& PC_idx){
    ifstream file("parameter/wp_setup/weld_top.txt");
    string data;
    mat row_tmp;
    mat PC_idx_tmp;
    int cnt = 0;
    while (getline(file, data)) {
        if (data.length() == 0)
            continue;
        setMatrix_arma(row_tmp, data);
        if (cnt == 0){
            PC_idx_tmp = row_tmp;
            cnt++;
        } else {
            PC_idx_tmp = join_cols(PC_idx_tmp, row_tmp);
            cnt++;
        }
    }
    PC_idx = PC_idx_tmp;
}

void loadWeldBottom(MatrixXd& weldbottom){
    fstream file("parameter/wp_setup/weld_bottom.txt");
    string data;
    inputfield line;
    map<string, vector<float> > dicts;
    MatrixXd row_tmp;
    MatrixXd weld_bottom_tmp;
    int cnt =0;
    while (std::getline(file,data)) {
        // skip the empty line 
        if (data.length()==0)
            continue;
        // customized matrix setting
        // get the matrix elements in string
        setMatrix(row_tmp, data);
        // concantenate current row with old weld traj
        if (cnt == 0){
            weld_bottom_tmp = row_tmp;
            cnt++;
        } // the first row 
        else {
            weld_bottom_tmp = Vcat(weld_bottom_tmp, row_tmp);
            cnt++;
        }
    }
    weldbottom = weld_bottom_tmp;
}

void loadWeldBottom_arma(mat& PC_idx){
    ifstream file("parameter/wp_setup/weld_bottom.txt");
    string data;
    mat row_tmp;
    mat PC_idx_tmp;
    int cnt = 0;
    while (getline(file, data)) {
        if (data.length() == 0)
            continue;
        setMatrix_arma(row_tmp, data);
        if (cnt == 0){
            PC_idx_tmp = row_tmp;
            cnt++;
        } else {
            PC_idx_tmp = join_cols(PC_idx_tmp, row_tmp);
            cnt++;
        }
    }
    PC_idx = PC_idx_tmp;
}


void loadWeldIn(MatrixXd& weldin){
    fstream file("parameter/wp_setup/weld_in.txt");
    string data;
    inputfield line;
    map<string, vector<float> > dicts;
    MatrixXd row_tmp;
    MatrixXd weld_in_tmp;
    int cnt =0;
    while (std::getline(file,data)) {
        // skip the empty line 
        if (data.length()==0)
            continue;
        // customized matrix setting
        // get the matrix elements in string
        setMatrix(row_tmp, data);
        // concantenate current row with old weld traj
        if (cnt == 0){
            weld_in_tmp = row_tmp;
            cnt++;
        } // the first row 
        else {
            weld_in_tmp = Vcat(weld_in_tmp, row_tmp);
            cnt++;
        }
    }
    weldin = weld_in_tmp;
}

void loadWeldIn_arma(mat& PC_idx){
    ifstream file("parameter/wp_setup/weld_in.txt");
    string data;
    mat row_tmp;
    mat PC_idx_tmp;
    int cnt = 0;
    while (getline(file, data)) {
        if (data.length() == 0)
            continue;
        setMatrix_arma(row_tmp, data);
        if (cnt == 0){
            PC_idx_tmp = row_tmp;
            cnt++;
        } else {
            PC_idx_tmp = join_cols(PC_idx_tmp, row_tmp);
            cnt++;
        }
    }
    PC_idx = PC_idx_tmp;
}

void load_PC_idx(MatrixXd& PC_idx){
    fstream file("parameter/PC_idx.txt");
    string data;
    inputfield line;
    map<string, vector<float> > dicts;
    MatrixXd row_tmp;
    MatrixXd PC_idx_tmp;
    int cnt =0;
    while (std::getline(file,data)) {
        // skip the empty line 
        if (data.length()==0)
            continue;
        // customized matrix setting
        // get the matrix elements in string
        setMatrix(row_tmp, data);
        // concantenate current row with old weld traj
        if (cnt == 0){
            PC_idx_tmp = row_tmp;
            cnt++;
        } // the first row 
        else {
            PC_idx_tmp = Vcat(PC_idx_tmp, row_tmp);
            cnt++;
        }
    }
    PC_idx = PC_idx_tmp;
}

void load_PC_idx_arma(mat& PC_idx){
    ifstream file("parameter/PC_idx.txt");
    string data;
    mat row_tmp;
    mat PC_idx_tmp;
    int cnt = 0;
    while (getline(file, data)) {
        if (data.length() == 0)
            continue;
        setMatrix_arma(row_tmp, data);
        if (cnt == 0){
            PC_idx_tmp = row_tmp;
            cnt++;
        } else {
            PC_idx_tmp = join_cols(PC_idx_tmp, row_tmp);
            cnt++;
        }
    }
    PC_idx = PC_idx_tmp;
}


void loadWeldOut(MatrixXd& weldout){
    fstream file("parameter/wp_setup/weld_out.txt");
    string data;
    inputfield line;
    map<string, vector<float> > dicts;
    MatrixXd row_tmp;
    MatrixXd weld_out_tmp;
    int cnt =0;
    while (std::getline(file,data)) {
        // skip the empty line 
        if (data.length()==0)
            continue;
        // customized matrix setting
        // get the matrix elements in string
        setMatrix(row_tmp, data);
        // concantenate current row with old weld traj
        if (cnt == 0){
            weld_out_tmp = row_tmp;
            cnt++;
        } // the first row 
        else {
            weld_out_tmp = Vcat(weld_out_tmp, row_tmp);
            cnt++;
        }
    }
    weldout = weld_out_tmp;
}


void loadWeldOut_arma(mat& PC_idx){
    ifstream file("parameter/wp_setup/weld_out.txt");
    string data;
    mat row_tmp;
    mat PC_idx_tmp;
    int cnt = 0;
    while (getline(file, data)) {
        if (data.length() == 0)
            continue;
        setMatrix_arma(row_tmp, data);
        if (cnt == 0){
            PC_idx_tmp = row_tmp;
            cnt++;
        } else {
            PC_idx_tmp = join_cols(PC_idx_tmp, row_tmp);
            cnt++;
        }
    }
    PC_idx = PC_idx_tmp;
}


void loadplanePoints(MatrixXd& planePoints){
    ifstream file("parameter/wp_setup/planePoints.txt");
    string data;
    inputfield line;
    map<string, vector<float> > dicts;
    MatrixXd row_tmp;
    MatrixXd pp_tmp;
    int cnt = 0;
    while (std::getline(file, data)) {
        // skip the empty line 
        if (data.length()==0)
            continue;
        // customized matrix setting
        // get the matrix elements in string
        setMatrix(row_tmp, data);
        // concantenate current row with old weld traj
        if (cnt == 0){
            pp_tmp = row_tmp;
            cnt++;
        } // the first row  
        else{
            pp_tmp = Vcat(pp_tmp, row_tmp);
            cnt++;
        }  
    }
    // set the original weldTraj
    planePoints = pp_tmp;
}

void loadplanePoints_arma(mat& PC_idx){
    ifstream file("parameter/wp_setup/planePoints.txt");
    string data;
    mat row_tmp;
    mat PC_idx_tmp;
    int cnt = 0;
    while (getline(file, data)) {
        if (data.length() == 0)
            continue;
        setMatrix_arma(row_tmp, data);
        if (cnt == 0){
            PC_idx_tmp = row_tmp;
            cnt++;
        } else {
            PC_idx_tmp = join_cols(PC_idx_tmp, row_tmp);
            cnt++;
        }
    }
    PC_idx = PC_idx_tmp;
}

void loadabc(MatrixXd& abc){
    ifstream file("parameter/wp_setup/abc.txt");
    string data;
    inputfield line;
    map<string, vector<float> > dicts;
    MatrixXd row_tmp;
    MatrixXd abc_tmp;
    int cnt = 0;
    while (std::getline(file, data)) {
        // skip the empty line 
        if (data.length()==0)
            continue;
        // customized matrix setting
        // get the matrix elements in string
        setMatrix(row_tmp, data);
        // concantenate current row with old weld traj
        if (cnt == 0){
            abc_tmp = row_tmp;
            cnt++;
        } // the first row  
        else{
            abc_tmp = Vcat(abc_tmp, row_tmp);
            cnt++;
        }  
    }
    // set the original weldTraj
    abc = abc_tmp;
}

void loadabc_arma(mat& PC_idx){
    ifstream file("parameter/wp_setup/abc.txt");
    string data;
    mat row_tmp;
    mat PC_idx_tmp;
    int cnt = 0;
    while (getline(file, data)) {
        if (data.length() == 0)
            continue;
        setMatrix_arma(row_tmp, data);
        if (cnt == 0){
            PC_idx_tmp = row_tmp;
            cnt++;
        } else {
            PC_idx_tmp = join_cols(PC_idx_tmp, row_tmp);
            cnt++;
        }
    }
    PC_idx = PC_idx_tmp;
}

void loadArr(MatrixXd& weldTraj){
    /*
    * load the weld points trajectories 
    * weldTraj is stored as double type matrix 
    * size: 3*n (n is the weld point number)
    * input format: n*3 (n is weld point number, 3 is x,y,z in sequence)
    */
    ifstream file("parameter/wp_setup/real_weld.txt");
    string data;
    inputfield line;
    map<string, vector<float> > dicts;
    MatrixXd row_tmp;
    MatrixXd weldTraj_tmp;
    int cnt = 0;
    cout << "start loading weld trajctory points:" << endl;
    cout << "..." << endl;
    while (std::getline(file, data)) {
        // skip the empty line 
        if (data.length()==0)
            continue;
        // customized matrix setting
        // get the matrix elements in string
        setMatrix(row_tmp, data);
        // concantenate current row with old weld traj
        if (cnt == 0){
            weldTraj_tmp = row_tmp;
            cnt++;
        } // the first row  
        else{
            weldTraj_tmp = Vcat(weldTraj_tmp, row_tmp);
            cnt++;
        }  
    }
    // set the original weldTraj
    weldTraj = weldTraj_tmp;
    // cout << "loaded weldTraj points are:" << endl;
    // cout << weldTraj << endl;
}


void loadAnchor(MatrixXd& point_anchor){
    ifstream file("parameter/wp_setup/point_anchor.txt");
    string data;
    inputfield line;
    map<string, vector<float> > dicts;
    MatrixXd row_tmp;
    MatrixXd pa_tmp;
    int cnt = 0;
    while (std::getline(file, data)) {
        // skip the empty line 
        if (data.length()==0)
            continue;
        // customized matrix setting
        // get the matrix elements in string
        setMatrix(row_tmp, data);
        // concantenate current row with old weld traj
        if (cnt == 0){
            pa_tmp = row_tmp;
            cnt++;
        } // the first row  
        else{
            pa_tmp = Vcat(pa_tmp, row_tmp);
            cnt++;
        }  
    }
    // set the original weldTraj
    point_anchor = pa_tmp;
}

void loadAnchor_arma(mat& PC_idx){
    ifstream file("parameter/wp_setup/point_anchor.txt");
    string data;
    mat row_tmp;
    mat PC_idx_tmp;
    int cnt = 0;
    while (getline(file, data)) {
        if (data.length() == 0)
            continue;
        setMatrix_arma(row_tmp, data);
        if (cnt == 0){
            PC_idx_tmp = row_tmp;
            cnt++;
        } else {
            PC_idx_tmp = join_cols(PC_idx_tmp, row_tmp);
            cnt++;
        }
    }
    PC_idx = PC_idx_tmp;
}

void loadDHbase(MatrixXd &DH, MatrixXd &base){
    YAML::Node config = YAML::LoadFile("parameter/Kinematics.yaml");
    //cout<< "Test:" << config["DH_params"].Type()<<endl;
    //cout<< "Test:" << config.Type()<<endl;
    assert(config["DH_params"].Type() == YAML::NodeType::Sequence);
    assert(config.Type() == YAML::NodeType::Map);

    // discard the first four DH parameters, we will use it as base 
    // use the remaining DH parameters 
    int cout = 1;
    std::vector<double> remain_DH;
    double tmp;

    // initialize base 
    base.resize(3,1); // base is x,y,z : 3-dimensional vector
    for (YAML::const_iterator it=config["DH_params"].begin(); it!=config["DH_params"].end(); ++it){
        // the dependes on the DH parameters property, z offset is d entry of DH parameter 
        // first dimension 
        if (cout == 2){
            base(0,0) = 0;
            base(1,0) = 0;
            base(2,0) = it->as<double>();
        }
        if (cout >= 5){
            // tmp = it->as<double>();
            // std::cout << tmp << std::endl;
            // abort();
            remain_DH.push_back(it->as<double>());
        }
        cout++;
    }
    int row, col;
    col = 4; // this is default format for DH 
    row = remain_DH.size()/col;
    vector2Matrix(DH, row, col, remain_DH);
}

void loadDHbase_arma(arma::mat& DH, arma::mat& base) {
    YAML::Node config = YAML::LoadFile("parameter/Kinematics.yaml");
    assert(config["DH_params"].Type() == YAML::NodeType::Sequence);
    assert(config.Type() == YAML::NodeType::Map);

    int cout = 1;
    std::vector<double> remain_DH;

    // Initialize base
    base.resize(3, 1);
    base.fill(0.0); // Initialize with zeros
    for (YAML::const_iterator it = config["DH_params"].begin(); it != config["DH_params"].end(); ++it) {
        if (cout == 2) {
            base(2, 0) = it->as<double>(); // Set the z offset
        }
        if (cout >= 5) {
            remain_DH.push_back(it->as<double>());
        }
        cout++;
    }

    int row, col;
    col = 4; // This is the default format for DH parameters
    row = remain_DH.size() / col;

    // Convert vector to Armadillo matrix
    DH = arma::mat(row, col);
    for (int r = 0; r < row; ++r) {
        for (int c = 0; c < col; ++c) {
            DH(r, c) = remain_DH[r * col + c];
        }
    }
}

void setInputField(inputfield& line, string& data){
    /*
    * make sure each txt line has no ' ' in the end
    * make sure there is '\n' at end of txt file 
    */
    // words count 
    int cnt = 0;
    string word = ""; 
    for (auto x : data){
        if (x == ' ') 
       {    
            if(cnt == 0){
                // get field
                line.field.assign(word);
                word = ""; // empty word
                cnt++;
            }
            else{
                // get elements
                // cout << word << endl;
                line.elements.push_back(stof(word)); // push back next element in float
                word = ""; // empty word
                cnt++;
            }
       } 
       else
       { 
           word = word + x; 
       } 
    }
    // the last element doesn't have blanket behind it
    // cout << word << endl;
    line.elements.push_back(stof(word)); // push back next element in float

    // print out to test the correctness
    // cout << line.field << endl;
    // for (vector<float>::const_iterator i = line.elements.begin(); i != line.elements.end(); ++i)
    //     cout << *i << ' ';
    // cout << endl <<  "good to go" << endl;
}

struct Point{
float x, y, z;
};


void loadWeldTraj(MatrixXd& weldTraj){
    /*
    * load the weld points trajectories 
    * weldTraj is stored as double type matrix 
    * size: 3*n (n is the weld point number)
    * input format: n*3 (n is weld point number, 3 is x,y,z in sequence)
    */
    // ifstream file("src/pointcloud/weldTraj.txt");
    fstream file("parameter/wp_setup/weldTraj.txt");
    string data;
    inputfield line;
    map<string, vector<float> > dicts;
    MatrixXd row_tmp;
    MatrixXd weldTraj_tmp;
    int cnt = 0;
    while (std::getline(file, data)) {
        // skip the empty line 
        if (data.length()==0)
            continue;
        // customized matrix setting
        // get the matrix elements in string
        setMatrix(row_tmp, data);
        // concantenate current row with old weld traj
        if (cnt == 0){
            weldTraj_tmp = row_tmp;
            cnt++;
        } // the first row  
        else{
            weldTraj_tmp = Vcat(weldTraj_tmp, row_tmp);
            cnt++;
        }  
    }
    // set the original weldTraj
    weldTraj = weldTraj_tmp;
    // cout << "loaded weldTraj points are:" << endl;
    // cout << weldTraj << endl;
}

void loadWeldTraj_arma(mat& weldTraj){
    /*
    * load the weld points trajectories 
    * weldTraj is stored as double type matrix 
    * size: 3*n (n is the weld point number)
    * input format: n*3 (n is weld point number, 3 is x,y,z in sequence)
    */
    // ifstream file("src/pointcloud/weldTraj.txt");
    fstream file("parameter/wp_setup/weldTraj.txt");
    string data;
    mat row_tmp;
    mat weldTraj_tmp;
    int cnt = 0;
    while (std::getline(file, data)) {
        // skip the empty line 
        if (data.length()==0)
            continue;
        // customized matrix setting
        // get the matrix elements in string
        setMatrix_arma(row_tmp, data);
        // concantenate current row with old weld traj
        if (cnt == 0){
            weldTraj_tmp = row_tmp;
            cnt++;
        } // the first row  
        else{
            weldTraj_tmp = join_cols(weldTraj_tmp, row_tmp);
            cnt++;
        }  
    }
    // set the original weldTraj
    weldTraj = weldTraj_tmp;
    // cout << "loaded weldTraj points are:" << endl;
    // cout << weldTraj << endl;
}

void loadWorkpieceSetting(MatrixXd& abc,  MatrixXd& planePoints, MatrixXd& point_anchor, 
        MatrixXd& weld_bottom, MatrixXd& weld_in, MatrixXd& weld_left, MatrixXd& weld_out, MatrixXd& weld_right, MatrixXd& weld_top){
    loadabc(abc);
    loadAnchor(point_anchor);
    loadplanePoints(planePoints);
    loadWeldIn(weld_in);
    loadWeldBottom(weld_bottom);
    loadWeldLeft(weld_left);
    loadWeldOut(weld_out);
    loadWeldRight(weld_right);
    loadWeldTop(weld_top);
}

void loadWorkpieceSetting_arma(mat& abc,  mat& planePoints, mat& point_anchor, 
        mat& weld_bottom, mat& weld_in, mat& weld_left, mat& weld_out, mat& weld_right, mat& weld_top){
    loadabc_arma(abc);
    loadAnchor_arma(point_anchor);
    loadplanePoints_arma(planePoints);
    loadWeldIn_arma(weld_in);
    loadWeldBottom_arma(weld_bottom);
    loadWeldLeft_arma(weld_left);
    loadWeldOut_arma(weld_out);
    loadWeldRight_arma(weld_right);
    loadWeldTop_arma(weld_top);
}


int loadSafePolishingSetting(double& alphaY_limit, double& alphaZ_limit, double& thres, int& nstep1, int& nstep2, int& nwait, MatrixXd& theta_init_default, MatrixXd& wp_pos_init_default){
    ifstream file("parameter/parameters.txt");
    string data;
    inputfield line;
    map<string, vector<float> > dicts;
    while (std::getline(file, data)) {
        // cout << data << "\n";
        // clear line 
        line.field.clear();
        line.elements.clear();
        
        // customized parameter setting
        // get the elements in string
        // the first elements is name field
        setInputField(line, data);

        // feed line into dictionary 
        dicts.insert(pair<string, vector<float> >(line.field, line.elements));
    }
    // set the input value
    alphaY_limit = (double) dicts["alphaY_limit"][0];
    alphaZ_limit = (double) dicts["alphaZ_limit"][0];
    // theta = (double) dicts["theta"][0];
    // start2exe = (int) dicts["start2exe"][0];
    nstep1 = (int) dicts["nstep1"][0];
    nstep2 = (int) dicts["nstep2"][0];
    nwait = (int) dicts["nwait"][0];
    thres = (double) dicts["thres"][0];

    theta_init_default << (double) dicts["theta_init_default"][0],
                  (double) dicts["theta_init_default"][1],
                  (double) dicts["theta_init_default"][2],
                  (double) dicts["theta_init_default"][3],
                  (double) dicts["theta_init_default"][4],
                  (double) dicts["theta_init_default"][5];

    wp_pos_init_default << (double) dicts["wp_pos_init_default"][0],
                  (double) dicts["wp_pos_init_default"][1],
                  (double) dicts["wp_pos_init_default"][2],

    // check the input parameter setting value
    cout << "alphaY_limit is :" << alphaY_limit << endl;
    cout << "alphaZ_limit is :" << alphaZ_limit << endl;
    cout << "theta_init_default is :\n" << theta_init_default << endl;
    cout << "wp_pos_init_default is :\n" << wp_pos_init_default << endl;
    // set the input value
    cout << "all seems good" << endl;
    // abort();
    return 0;
}

int loadSafePolishingSetting_arma(double& alphaY_limit, double& alphaZ_limit, double& thres, int& nstep1, int& nstep2, int& nwait, mat& theta_init_default, mat& wp_pos_init_default){
    ifstream file("parameter/parameters.txt");
    string data;
    inputfield line;
    map<string, vector<float> > dicts;
    while (std::getline(file, data)) {
        // cout << data << "\n";
        // clear line 
        line.field.clear();
        line.elements.clear();
        
        // customized parameter setting
        // get the elements in string
        // the first elements is name field
        setInputField(line, data);

        // feed line into dictionary 
        dicts.insert(pair<string, vector<float> >(line.field, line.elements));
    }
    // set the input value
    alphaY_limit = (double) dicts["alphaY_limit"][0];
    alphaZ_limit = (double) dicts["alphaZ_limit"][0];
    // theta = (double) dicts["theta"][0];
    // start2exe = (int) dicts["start2exe"][0];
    nstep1 = (int) dicts["nstep1"][0];
    nstep2 = (int) dicts["nstep2"][0];
    nwait = (int) dicts["nwait"][0];
    thres = (double) dicts["thres"][0];

    theta_init_default << (double) dicts["theta_init_default"][0] << endr
                  << (double) dicts["theta_init_default"][1] << endr
                  << (double) dicts["theta_init_default"][2] << endr
                  << (double) dicts["theta_init_default"][3] << endr
                  << (double) dicts["theta_init_default"][4] << endr 
                  << (double) dicts["theta_init_default"][5] << endr;

    wp_pos_init_default << (double) dicts["wp_pos_init_default"][0] << endr
                  << (double) dicts["wp_pos_init_default"][1] << endr
                  << (double) dicts["wp_pos_init_default"][2] << endr;

    // check the input parameter setting value
    cout << "alphaY_limit is :" << alphaY_limit << endl;
    cout << "alphaZ_limit is :" << alphaZ_limit << endl;
    cout << "theta_init_default is :\n" << theta_init_default << endl;
    cout << "wp_pos_init_default is :\n" << wp_pos_init_default << endl;
    // set the input value
    cout << "all seems good" << endl;
    // abort();
    return 0;
}

void loadM1(MatrixXd& M){
    YAML::Node config = YAML::LoadFile("parameter/Kinematics.yaml");
    assert(config["M1"].Type() == YAML::NodeType::Sequence);
    assert(config.Type() == YAML::NodeType::Map);
    int cout = 1;
    std::vector<double> M1_temp;
    // record data to vector
    for (YAML::const_iterator it=config["M1"].begin(); it!=config["M1"].end(); ++it){
        // the dependes on the DH parameters property, z offset is d entry of DH parameter 
        // first dimension 
        M1_temp.push_back(it->as<double>());
        cout++;
    }
    int row, col;
    col = 4; // this is default format for transformation matrix 
    row = M1_temp.size()/col;
    vector2Matrix(M, row, col, M1_temp);
}

void loadM1_arma(arma::mat& M) {
    YAML::Node config = YAML::LoadFile("parameter/Kinematics.yaml");
    assert(config["M1"].Type() == YAML::NodeType::Sequence);
    assert(config.Type() == YAML::NodeType::Map);
    
    int cout = 1;
    std::vector<double> M1_temp;

    // Record data to vector
    for (YAML::const_iterator it = config["M1"].begin(); it != config["M1"].end(); ++it) {
        // The depends on the DH parameters property, z offset is d entry of DH parameter
        // First dimension
        M1_temp.push_back(it->as<double>());
        cout++;
    }

    int row, col;
    col = 4; // This is the default format for a transformation matrix
    row = M1_temp.size() / col;

    // Convert vector to Armadillo matrix
    M = arma::mat(row, col);
    for (int r = 0; r < row; ++r) {
        for (int c = 0; c < col; ++c) {
            M(r, c) = M1_temp[r * col + c];
        }
    }
}

void loadjnt2tool(MatrixXd& M){
    YAML::Node config = YAML::LoadFile("parameter/Kinematics.yaml");
    assert(config["J6_to_Tool_M"].Type() == YAML::NodeType::Sequence);
    assert(config.Type() == YAML::NodeType::Map);
    int cout = 1;
    std::vector<double> j2tool;
    // record data to vector
    for (YAML::const_iterator it=config["J6_to_Tool_M"].begin(); it!=config["J6_to_Tool_M"].end(); ++it){
        // the dependes on the DH parameters property, z offset is d entry of DH parameter 
        // first dimension 
        j2tool.push_back(it->as<double>());
        cout++;
    }
    int row, col;
    col = 4; // this is default format for transformation matrix 
    row = j2tool.size()/col;
    vector2Matrix(M, row, col, j2tool);
}

void loadjnt2tool_arma(arma::mat& M) {
    YAML::Node config = YAML::LoadFile("parameter/Kinematics.yaml");
    assert(config["J6_to_Tool_M"].Type() == YAML::NodeType::Sequence);
    assert(config.Type() == YAML::NodeType::Map);
    
    int cout = 1;
    std::vector<double> j2tool;

    // Record data to vector
    for (YAML::const_iterator it = config["J6_to_Tool_M"].begin(); it != config["J6_to_Tool_M"].end(); ++it) {
        // The depends on the DH parameters property, z offset is d entry of DH parameter
        // First dimension
        j2tool.push_back(it->as<double>());
        cout++;
    }

    int row, col;
    col = 4; // This is the default format for a transformation matrix
    row = j2tool.size() / col;

    // Convert vector to Armadillo matrix
    M = arma::mat(row, col);
    for (int r = 0; r < row; ++r) {
        for (int c = 0; c < col; ++c) {
            M(r, c) = j2tool[r * col + c];
        }
    }
}



void loadjnt2laser(MatrixXd& M){
    YAML::Node config = YAML::LoadFile("parameter/Kinematics.yaml");
    assert(config["J6_to_Laser_M"].Type() == YAML::NodeType::Sequence);
    assert(config.Type() == YAML::NodeType::Map);
    int cout = 1;
    std::vector<double> j2laser;
    // record data to vector
    for (YAML::const_iterator it=config["J6_to_Laser_M"].begin(); it!=config["J6_to_Laser_M"].end(); ++it){
        // the dependes on the DH parameters property, z offset is d entry of DH parameter 
        // first dimension 
        j2laser.push_back(it->as<double>());
        cout++;
    }
    int row, col;
    col = 4; // this is default format for transformation matrix 
    row = j2laser.size()/col;
    vector2Matrix(M, row, col, j2laser);
}



#endif