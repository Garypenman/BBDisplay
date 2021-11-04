#include "TString.h"
#include "TObjString.h"
#include "TObjArray.h"
#include <fstream>

#include "ConfigParser.h"

void ConfigParser(const char *configfilename){

  ifstream configfile(configfilename);
  

  if( configfile ){
    TString currentline;
    
    while( currentline.ReadLine(configfile) && !currentline.BeginsWith("endconfig")){
      if( !currentline.BeginsWith("#") ){
	TObjArray *tokens = currentline.Tokenize(" ");

	int ntokens = tokens->GetEntries();

	if( ntokens >= 2 ){
	  TString skey = ( (TObjString*) (*tokens)[0] )->GetString();

	  if( skey == "nlayers" ){
	    TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	    nlayers = stemp.Atoi();
	  }
	  
	  if( skey == "GEMtype" && ntokens >= nlayers + 1){
	    for(int i=1; i<ntokens; i++){
	      TString stemp = ( (TObjString*) (*tokens)[i] )->GetString();
	      GEMtype.push_back(stemp);
	    }
          }

	  if( skey == "nmodules" && ntokens >= nlayers + 1){
            for(int i=1; i<ntokens; i++){
              TString stemp = ( (TObjString*) (*tokens)[i] )->GetString();
              nmodules.push_back(stemp.Atoi());
            }
          }

	  if( skey == "U_strips" && ntokens >= nlayers + 1){
            for(int i=1; i<ntokens; i++){
              TString stemp = ( (TObjString*) (*tokens)[i] )->GetString();
              U_strips.push_back(stemp.Atof());
            }
          }

	  if( skey == "V_strips" && ntokens >= nlayers + 1){
            for(int i=1; i<ntokens; i++){
              TString stemp = ( (TObjString*) (*tokens)[i] )->GetString();
              V_strips.push_back(stemp.Atof());
            }
          }



	}
      }
    }
  }

  for(int i=0; i<nlayers; i++){
    if(GEMtype[i] == "none") continue;

    layer_list.push_back(i);
  }


}
