//TESTING THE FUNCTION read_weights

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <iostream>
using namespace std;


//particle
// 0 = Bu
// 1 = Bs

#define particle 0

//systematics
double read_weights(TString variable, double var_value);

double getWeight(double var_value, TH1D* h_weight);

int main(){

  double resultado;

  //read_weights(variable, var_value);

  for(int i=0; i<100; i++){
    resultado = read_weights("Bpt",i);
    cout<<i;
    cout << "weight" << resultado << endl;
  }

}
//main function ends

//definição da função:
double read_weights(TString variable, double var_value){
  
  TString input_file = particle ? "./results/Bs/mc_validation_plots/weights/weights.root" :  "./results/Bu/mc_validation_plots/weights/weights.root";

  TFile* f_wei = new TFile(input_file, "read");

  TH1D* histo_variable = (TH1D*)f_wei->Get(Form("weights_"+variable));

  double weight;
  double variable_min;
  double variable_max;

  variable_min = histo_variable->GetXaxis()->GetXmin();
  variable_max = histo_variable->GetXaxis()->GetXmax();  

  cout<<"min:"<<variable_min<<endl;
  cout<<"max:"<<variable_max<<endl;
  
  //if the event is not in the range its weight is 1.
  
  if(var_value>=variable_min && var_value<=variable_max){  
    weight = getWeight(var_value,histo_variable);
}
  else{
    weight = 1;
}


  return weight;
}

double getWeight(double var_value, TH1D* h_weight){
  int bin = h_weight->FindBin(var_value);
  return h_weight->GetBinContent(bin);
}
