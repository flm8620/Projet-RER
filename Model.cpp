#include <Imagine/Images.h>
#include <Imagine/Graphics.h>
#include <cmath>
#include "normaldistribution.h"
using namespace Imagine;
using namespace std;


const int nbportes = 30;
const int nbprevoyants = 1000;
const int nbconfort = 500;
const int tempsstation = 30;
const double propstresse = 0.8;
const int debitentree = 50;
const double ecarttype = 1;
const int nbdestinations = 3;
const double sigma = 5.0;
//tableau contenant le nombre de passagers dans le train
double s[30];
double T[300];
int Destination[nbdestinations];
double proportiondestination[nbdestinations];
const int nbsorties=2;
int sorties[nbsorties];
//proportion des passagers qui arrive par chaque sortie
double proportionsorties[nbsorties];
//paramètre de la loi exponentielle
const int lambda = -1;

void init_sorties(){
	/*int nbsorties;
	cout << "Donnez le nombre de sortie de la gare (entier entre 1 et 30)" << endl;
	cin >> nbsorties;
	int sorties[nbsorties];*/
	cout << "Designez la porte la plus pres de la sortie de la gare" << endl;
	cout << "Donnez la proportion de gens qui arrive par la sortie" << endl;
	double S=0;
	for (int i = 0; i < nbsorties; i++){
		cout << "pour la sortie" << i+1 << " (entier entre 1 et 30)" << endl;
		cin >> sorties[i];
		cout << "un double quelconque" << endl;
		cin >> proportionsorties[i];
		S += proportionsorties[i];
	}
	for (int i = 0; i < nbsorties; i++){
		proportionsorties[i] = proportionsorties[i] / S;
	}
}
void init_Destination(){
	cout << "Designez la porte la plus pres de chaque destination favorite" << endl;
	cout << "Donnez la proportion de gens qui veulent y aller" << endl;
	for (int i = 0; i < nbdestinations; i++){
		cout << "pour la porte" << i + 1 << " (un entier entre 1 et 30)" << endl;
		cin >> Destination[i];
		cout << "un double quelconque" << endl;
		cin >> proportiondestination[i];
	}
}

void gotominlocal(){
	for (int i = 0; i < nbsorties; i++){
		for (int j = 0; j < proportionsorties[i] * nbconfort; j++){
			if (sorties[i] == 1){
				for (int k = 1; k < nbportes; k++){
					if (s[k + 1] >= s[k]){ s[k] += 1; }
				}
			}
			if (sorties[i] == nbportes){
				for (int k = nbportes-1; k > 0; k--){
					if (s[k - 1] >= s[k]){ s[k] += 1; }
				}
			}
			if (sorties[i] < nbportes && sorties[i]>0){
                for (int k = sorties[i]; k < max( nbportes - sorties[i], sorties[i] ); k++){
					if (s[sorties[i] + k + 1] >= s[sorties[i] + k]){ s[sorties[i] + k] += 1; }
					else if (s[sorties[i] - k - 1] >= s[sorties[i] - k]){ s[sorties[i] - k] += 1; }
				}
			}
		}
	}
}

void loiexponentiel(int sorties[nbsorties]){
	//boucle sur toutes les sorties (1 loi par sortie)
	for (int i = 0; i < nbsorties; i++){
		//Cas sortie au début du quai
		if (sorties[i] == 0){
			for (int j = 0; j < nbportes; j++){
				s[j] += (1 - exp(lambda)) / (1 - exp(lambda*nbportes))*exp(lambda*j)*tempsstation*debitentree*propstresse*proportionsorties[i];
				//la loi est normée puis elle est multpliée par le nombre de passager stressé et en retard entrant par la sortie i.
			}
		}
		//Cas sortie à la fin du quai
		if (sorties[i] == nbportes - 1){
			for (int j = 0; j < nbportes; j++){
				s[j] += (1 - exp(lambda)) / (1 - exp(lambda*nbportes))*exp(lambda*(nbportes - j))*tempsstation*debitentree*propstresse*proportionsorties[i];
			}
		}
		//Cas sortie au milieu du quai
		if (sorties[i]>0 && sorties[i] < nbportes){
			for (int j = 0; j < sorties[i]; j++){
				s[j] += (1 - exp(lambda)) / (1 + exp(1) - exp(lambda*sorties[i]) - exp(lambda*(nbportes - sorties[i])))*exp(lambda*(sorties[i] - 1 - j))*tempsstation*debitentree*propstresse*proportionsorties[i];
			}
			for (int j = sorties[i]; j < nbportes; j++){
				s[j] += (1 - exp(lambda)) / (2 - exp(lambda*sorties[i]) - exp(lambda*(nbportes - sorties[i])))*exp(lambda*(j + 1 - sorties[i]))*tempsstation*debitentree*propstresse*proportionsorties[i];
			}
			//cette expression est la synthèse de trois étapes
			//1/ à la porte la plus proche de la sortie, la fonction a pour valeur 1 et elle décroit de chaque côté selon une même loi exponentielle de paramètre lambda
			//2/ je norme cette fonction sur le quai
			//3/ je multiplie par le nombre de passager stressé en retard entrant par la sortie i.
		}
	}
}

void loiuniforme(int sorties[nbsorties]){
	//boucle sur toutes les sorties
	for (int i = 0; i < nbsorties; i++){
		//boucle sur les destinations
		for (int j = 0; j < nbdestinations; j++){
			//boucle sur les portes entre la sortie et les destinations
			for (int k = min(i, j); k < max(i, j); k++){
				s[k] = tempsstation*debitentree*(1 - propstresse)*proportionsorties[i] * proportiondestination[j]/(j-i);
			}
		}
	}
}

void loinormale(int i, int j){
    //i - No.destination
    //j - No.porte
    NormalDistribution normal;
    double xSup=(j+1-Destination[i])/sigma;
    double xInf=(j - Destination[i])/sigma;
    double value=nbprevoyants*proportiondestination[i] * (normal.getPhi(xSup) - normal.getPhi(xInf));
    s[j] += value;
}

void init_loinormale(){
	for (int i = 0; i < 300; i++){ T[i] = 0; }
	T[0] = 0.5000; T[1] = 0.5040; T[2] = 0.5080; T[3] = 0.5120; T[4] = 0.5160; T[5] = 0.5199; T[6] = 0.5239; T[7] = 0.5279; T[8] = 0.5319; T[9] = 0.5359;
	T[10] = 0.5398; T[11] = 0.5438; T[12] = 0.5478; T[13] = 0.5517; T[14] = 0.5557; T[15] = 0.5596; T[16] = 0.5636; T[17] = 0.5675; T[18] = 0.5714; T[19] = 0.5753;
	T[20] = 0.5793; T[21] = 0.5832; T[22] = 0.5871; T[23] = 0.5910; T[24] = 0.5948; T[25] = 0.5987; T[26] = 0.6026; T[27] = 0.6064; T[28] = 0.6103; T[29] = 0.6141;
	T[30] = 0.6179; T[31] = 0.6217; T[32] = 0.6255; T[33] = 0.6293; T[34] = 0.6331; T[35] = 0.6368; T[36] = 0.6406; T[37] = 0.6443; T[38] = 0.6480; T[39] = 0.6517;
	T[40] = 0.6554; T[41] = 0.6591; T[42] = 0.6628; T[43] = 0.6664; T[44] = 0.6700; T[45] = 0.6736; T[46] = 0.6772; T[47] = 0.6808; T[48] = 0.6844; T[49] = 0.6879;
	T[50] = 0.6915; T[51] = 0.6950; T[52] = 0.6985; T[53] = 0.7019; T[54] = 0.7054; T[55] = 0.7088; T[56] = 0.7123; T[57] = 0.7157; T[58] = 0.7190; T[59] = 0.7224;
	T[60] = 0.7257; T[61] = 0.7291; T[62] = 0.7324; T[63] = 0.7357; T[64] = 0.7389; T[65] = 0.7422; T[66] = 0.7454; T[67] = 0.7486; T[68] = 0.7517; T[69] = 0.7549;
	T[70] = 0.7580; T[71] = 0.7611; T[72] = 0.7642; T[73] = 0.7673; T[74] = 0.7704; T[75] = 0.7734; T[76] = 0.7764; T[77] = 0.7794; T[78] = 0.7823; T[79] = 0.7852;
	T[80] = 0.7881; T[81] = 0.7910; T[82] = 0.7939; T[83] = 0.7967; T[84] = 0.7995; T[85] = 0.8023; T[86] = 0.8051; T[87] = 0.8078; T[88] = 0.8106; T[89] = 0.8133;
	T[90] = 0.8159; T[91] = 0.8186; T[92] = 0.8212; T[93] = 0.8238; T[94] = 0.8264; T[95] = 0.8289; T[96] = 0.8315; T[97] = 0.8340; T[98] = 0.8365; T[99] = 0.8389;
	T[100] = 0.8413; T[101] = 0.8438; T[102] = 0.8461; T[103] = 0.8485; T[104] = 0.8508; T[105] = 0.8531; T[106] = 0.8554; T[107] = 0.8577; T[108] = 0.8599; T[109] = 0.8621;
	T[110] = 0.8643; T[111] = 0.8665; T[112] = 0.8686; T[113] = 0.8708; T[114] = 0.8729; T[115] = 0.8749; T[116] = 0.8770; T[117] = 0.8790; T[118] = 0.8810; T[119] = 0.8830;
	T[120] = 0.8849; T[121] = 0.8869; T[122] = 0.8888; T[123] = 0.8907; T[124] = 0.8925; T[125] = 0.8944; T[126] = 0.8962; T[127] = 0.8980; T[128] = 0.8997; T[129] = 0.9015;
	T[130] = 0.9032; T[131] = 0.9049; T[132] = 0.9066; T[133] = 0.9082; T[134] = 0.9099; T[135] = 0.9115; T[136] = 0.9131; T[137] = 0.9147; T[138] = 0.9162; T[139] = 0.9177;
	T[140] = 0.9192; T[141] = 0.9207; T[142] = 0.9222; T[143] = 0.9236; T[144] = 0.9251; T[145] = 0.9265; T[146] = 0.9279; T[147] = 0.9292; T[148] = 0.9306; T[149] = 0.9319;
	T[150] = 0.9332; T[151] = 0.9345; T[152] = 0.9357; T[153] = 0.9370; T[154] = 0.9382; T[155] = 0.9394; T[156] = 0.9406; T[157] = 0.9418; T[158] = 0.9429; T[159] = 0.9441;
	T[160] = 0.9452; T[161] = 0.9463; T[162] = 0.9474; T[163] = 0.9484; T[164] = 0.9495; T[165] = 0.9505; T[166] = 0.9515; T[167] = 0.9525; T[168] = 0.9535; T[169] = 0.9545;
	T[170] = 0.9554; T[171] = 0.9564; T[172] = 0.9573; T[173] = 0.9582; T[174] = 0.9591; T[175] = 0.9599; T[176] = 0.9608; T[177] = 0.9616; T[178] = 0.9625; T[179] = 0.9633;
	T[180] = 0.9641; T[181] = 0.9649; T[182] = 0.9656; T[183] = 0.9664; T[184] = 0.9671; T[185] = 0.9678; T[186] = 0.9686; T[187] = 0.9693; T[188] = 0.9699; T[189] = 0.9706;
	T[190] = 0.9713; T[191] = 0.9719; T[192] = 0.9726; T[193] = 0.9732; T[194] = 0.9738; T[195] = 0.9744; T[196] = 0.9750; T[197] = 0.9756; T[198] = 0.9761; T[199] = 0.9767;
	T[200] = 0.9772; T[201] = 0.9778; T[202] = 0.9783; T[203] = 0.9788; T[204] = 0.9793; T[205] = 0.9798; T[206] = 0.9803; T[207] = 0.9808; T[208] = 0.9812; T[209] = 0.9817;
	T[210] = 0.9826; T[211] = 0.9826; T[212] = 0.9830; T[213] = 0.9834; T[214] = 0.9838; T[215] = 0.9842; T[216] = 0.9846; T[217] = 0.9850; T[218] = 0.9854; T[219] = 0.9857;
	T[220] = 0.9864; T[221] = 0.9864; T[222] = 0.9868; T[223] = 0.9871; T[224] = 0.9875; T[225] = 0.9878; T[226] = 0.9881; T[227] = 0.9884; T[228] = 0.9887; T[229] = 0.9890;
	T[230] = 0.9896; T[231] = 0.9896; T[232] = 0.9898; T[233] = 0.9901; T[234] = 0.9904; T[235] = 0.9906; T[236] = 0.9909; T[237] = 0.9911; T[238] = 0.9913; T[239] = 0.9916;
	T[240] = 0.9920; T[241] = 0.9920; T[242] = 0.9922; T[243] = 0.9925; T[244] = 0.9927; T[245] = 0.9929; T[246] = 0.9931; T[247] = 0.9932; T[248] = 0.9934; T[249] = 0.9936;
	T[250] = 0.9940; T[251] = 0.9940; T[252] = 0.9941; T[253] = 0.9943; T[254] = 0.9945; T[255] = 0.9946; T[256] = 0.9948; T[257] = 0.9949; T[258] = 0.9951; T[259] = 0.9952;
	T[260] = 0.9955; T[261] = 0.9955; T[262] = 0.9956; T[263] = 0.9957; T[264] = 0.9959; T[265] = 0.9960; T[266] = 0.9961; T[267] = 0.9962; T[268] = 0.9963; T[269] = 0.9964;
	T[270] = 0.9966; T[271] = 0.9966; T[272] = 0.9967; T[273] = 0.9968; T[274] = 0.9969; T[275] = 0.9970; T[276] = 0.9971; T[277] = 0.9972; T[278] = 0.9973; T[279] = 0.9974;
	T[280] = 0.9975; T[281] = 0.9975; T[282] = 0.9976; T[283] = 0.9977; T[284] = 0.9977; T[285] = 0.9978; T[286] = 0.9979; T[287] = 0.9979; T[288] = 0.9980; T[289] = 0.9981;
	T[290] = 0.9982; T[291] = 0.9982; T[292] = 0.9982; T[293] = 0.9983; T[294] = 0.9984; T[295] = 0.9984; T[296] = 0.9985; T[297] = 0.9985; T[298] = 0.9986; T[299] = 1;

}


void init_quai(){
	for (int i = 0; i < nbportes; i++){
		s[i] = 0;
	}
}


int main(){
	init_quai();
	init_loinormale();
	init_Destination();
	init_sorties();

	for (int i = 0; i < nbdestinations; i++){
		for (int j = 0; j < nbportes; j++){
			loinormale(i, j);
		}
	}
	/**gotominlocal();
	for (int i = 0; i < nbsorties; i++){
		for (int j = 0; j < proportionsorties[i] * propstresse*tempsstation / debitentree; j++){
			s[sorties[i]] += 1;
		}
	}
	for (int i = 0; i < nbsorties; i++){
		for (int j = 0; j < nbdestinations; j++){
			for (int k = 0; k < abs(sorties[i] - Destination[j]); j++){
				s[k + sorties[i]] += proportionsorties[i] * (1 - )*proportiondestination[j] * tempsstation / debitentree;
			}
		}
	}**/
	loiexponentiel(sorties);
	loiuniforme(sorties);


	openWindow(600, 600);
	for (int i = 0; i < nbportes; i++){
		drawLine(20 * i, 600, 20 * i, 600.-s[i], MAGENTA);
		cout << s[i] << endl;
	}
	system("pause");
	return 0;
}
