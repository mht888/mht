#include "emp-tool/circuits/float32.h"
using emp::Float;
using emp::Bit;

Bit Float::less_than(const Float & rhs) const {
	Float res(*this);
	Bit *B = new Bit[378];

	B[0] = res[0];
	B[1] = res[1];
	B[2] = res[2];
	B[3] = res[3];
	B[4] = res[4];
	B[5] = res[5];
	B[6] = res[6];
	B[7] = res[7];
	B[8] = res[8];
	B[9] = res[9];
	B[10] = res[10];
	B[11] = res[11];
	B[12] = res[12];
	B[13] = res[13];
	B[14] = res[14];
	B[15] = res[15];
	B[16] = res[16];
	B[17] = res[17];
	B[18] = res[18];
	B[19] = res[19];
	B[20] = res[20];
	B[21] = res[21];
	B[22] = res[22];
	B[23] = res[23];
	B[24] = res[24];
	B[25] = res[25];
	B[26] = res[26];
	B[27] = res[27];
	B[28] = res[28];
	B[29] = res[29];
	B[30] = res[30];
	B[31] = res[31];
	B[32] = rhs[0];
	B[33] = rhs[1];
	B[34] = rhs[2];
	B[35] = rhs[3];
	B[36] = rhs[4];
	B[37] = rhs[5];
	B[38] = rhs[6];
	B[39] = rhs[7];
	B[40] = rhs[8];
	B[41] = rhs[9];
	B[42] = rhs[10];
	B[43] = rhs[11];
	B[44] = rhs[12];
	B[45] = rhs[13];
	B[46] = rhs[14];
	B[47] = rhs[15];
	B[48] = rhs[16];
	B[49] = rhs[17];
	B[50] = rhs[18];
	B[51] = rhs[19];
	B[52] = rhs[20];
	B[53] = rhs[21];
	B[54] = rhs[22];
	B[55] = rhs[23];
	B[56] = rhs[24];
	B[57] = rhs[25];
	B[58] = rhs[26];
	B[59] = rhs[27];
	B[60] = rhs[28];
	B[61] = rhs[29];
	B[62] = rhs[30];
	B[63] = rhs[31];
	B[306] = !B[12];
	B[64] = B[44] & B[306];
	B[64] = !B[64];
	B[307] = !B[20];
	B[65] = B[52] & B[307];
	B[65] = !B[65];
	B[308] = !B[8];
	B[66] = B[40] & B[308];
	B[66] = !B[66];
	B[309] = !B[4];
	B[67] = B[36] & B[309];
	B[67] = !B[67];
	B[310] = !B[16];
	B[68] = B[48] & B[310];
	B[68] = !B[68];
	B[311] = !B[1];
	B[69] = B[33] & B[311];
	B[69] = !B[69];
	B[70] = B[24] | B[23];
	B[70] = !B[70];
	B[71] = !B[26];
	B[72] = !B[25];
	B[73] = B[71] & B[72];
	B[74] = B[70] & B[73];
	B[75] = !B[28];
	B[76] = !B[27];
	B[77] = B[75] & B[76];
	B[78] = !B[30];
	B[79] = !B[29];
	B[80] = B[78] & B[79];
	B[81] = B[77] & B[80];
	B[82] = B[74] & B[81];
	B[312] = !B[82];
	B[83] = B[31] & B[312];
	B[84] = !B[56];
	B[85] = !B[55];
	B[86] = B[84] & B[85];
	B[87] = !B[58];
	B[88] = !B[57];
	B[89] = B[87] & B[88];
	B[90] = B[86] & B[89];
	B[91] = !B[60];
	B[92] = !B[59];
	B[93] = B[91] & B[92];
	B[94] = !B[62];
	B[95] = !B[61];
	B[96] = B[94] & B[95];
	B[97] = B[93] & B[96];
	B[98] = B[90] & B[97];
	B[98] = !B[98];
	B[99] = B[63] & B[98];
	B[100] = B[83] ^ B[99];
	B[100] = !B[100];
	B[313] = !B[21];
	B[101] = B[53] & B[313];
	B[101] = !B[101];
	B[314] = !B[52];
	B[102] = B[101] & B[314];
	B[103] = B[20] & B[102];
	B[103] = !B[103];
	B[315] = !B[53];
	B[104] = B[21] & B[315];
	B[104] = !B[104];
	B[105] = B[103] & B[104];
	B[316] = !B[49];
	B[106] = B[17] & B[316];
	B[106] = !B[106];
	B[317] = !B[17];
	B[107] = B[49] & B[317];
	B[107] = !B[107];
	B[318] = !B[48];
	B[108] = B[107] & B[318];
	B[109] = B[16] & B[108];
	B[109] = !B[109];
	B[110] = B[106] & B[109];
	B[110] = !B[110];
	B[319] = !B[18];
	B[111] = B[50] & B[319];
	B[111] = !B[111];
	B[320] = !B[19];
	B[112] = B[51] & B[320];
	B[112] = !B[112];
	B[113] = B[111] & B[112];
	B[113] = !B[113];
	B[321] = !B[113];
	B[114] = B[110] & B[321];
	B[114] = !B[114];
	B[322] = !B[51];
	B[115] = B[19] & B[322];
	B[115] = !B[115];
	B[116] = B[114] & B[115];
	B[323] = !B[50];
	B[117] = B[112] & B[323];
	B[118] = B[117] & B[18];
	B[118] = !B[118];
	B[119] = B[116] & B[118];
	B[119] = !B[119];
	B[120] = B[65] & B[101];
	B[120] = !B[120];
	B[324] = !B[120];
	B[121] = B[119] & B[324];
	B[121] = !B[121];
	B[122] = B[105] & B[121];
	B[122] = !B[122];
	B[325] = !B[22];
	B[123] = B[54] & B[325];
	B[326] = !B[123];
	B[124] = B[122] & B[326];
	B[124] = !B[124];
	B[327] = !B[54];
	B[125] = B[22] & B[327];
	B[125] = !B[125];
	B[126] = B[124] & B[125];
	B[328] = !B[47];
	B[127] = B[15] & B[328];
	B[127] = !B[127];
	B[329] = !B[45];
	B[128] = B[13] & B[329];
	B[128] = !B[128];
	B[330] = !B[13];
	B[129] = B[45] & B[330];
	B[129] = !B[129];
	B[331] = !B[44];
	B[130] = B[129] & B[331];
	B[131] = B[12] & B[130];
	B[131] = !B[131];
	B[132] = B[128] & B[131];
	B[132] = !B[132];
	B[332] = !B[14];
	B[133] = B[46] & B[332];
	B[133] = !B[133];
	B[333] = !B[15];
	B[134] = B[47] & B[333];
	B[134] = !B[134];
	B[135] = B[133] & B[134];
	B[135] = !B[135];
	B[334] = !B[135];
	B[136] = B[132] & B[334];
	B[136] = !B[136];
	B[137] = B[127] & B[136];
	B[335] = !B[46];
	B[138] = B[134] & B[335];
	B[139] = B[138] & B[14];
	B[139] = !B[139];
	B[336] = !B[41];
	B[140] = B[9] & B[336];
	B[140] = !B[140];
	B[337] = !B[9];
	B[141] = B[41] & B[337];
	B[141] = !B[141];
	B[338] = !B[40];
	B[142] = B[141] & B[338];
	B[143] = B[8] & B[142];
	B[143] = !B[143];
	B[144] = B[140] & B[143];
	B[144] = !B[144];
	B[339] = !B[10];
	B[145] = B[42] & B[339];
	B[145] = !B[145];
	B[340] = !B[11];
	B[146] = B[43] & B[340];
	B[146] = !B[146];
	B[147] = B[145] & B[146];
	B[147] = !B[147];
	B[341] = !B[147];
	B[148] = B[144] & B[341];
	B[148] = !B[148];
	B[342] = !B[43];
	B[149] = B[11] & B[342];
	B[149] = !B[149];
	B[150] = B[148] & B[149];
	B[343] = !B[42];
	B[151] = B[146] & B[343];
	B[152] = B[151] & B[10];
	B[152] = !B[152];
	B[153] = B[150] & B[152];
	B[153] = !B[153];
	B[154] = B[64] & B[129];
	B[344] = !B[135];
	B[155] = B[154] & B[344];
	B[155] = !B[155];
	B[345] = !B[155];
	B[156] = B[153] & B[345];
	B[156] = !B[156];
	B[157] = B[139] & B[156];
	B[158] = B[155] | B[147];
	B[158] = !B[158];
	B[159] = B[66] & B[141];
	B[346] = !B[37];
	B[160] = B[5] & B[346];
	B[160] = !B[160];
	B[347] = !B[5];
	B[161] = B[37] & B[347];
	B[161] = !B[161];
	B[348] = !B[36];
	B[162] = B[161] & B[348];
	B[163] = B[4] & B[162];
	B[163] = !B[163];
	B[164] = B[160] & B[163];
	B[164] = !B[164];
	B[349] = !B[6];
	B[165] = B[38] & B[349];
	B[165] = !B[165];
	B[350] = !B[7];
	B[166] = B[39] & B[350];
	B[166] = !B[166];
	B[167] = B[165] & B[166];
	B[167] = !B[167];
	B[351] = !B[167];
	B[168] = B[164] & B[351];
	B[168] = !B[168];
	B[352] = !B[39];
	B[169] = B[7] & B[352];
	B[169] = !B[169];
	B[170] = B[168] & B[169];
	B[353] = !B[38];
	B[171] = B[166] & B[353];
	B[172] = B[171] & B[6];
	B[172] = !B[172];
	B[354] = !B[2];
	B[173] = B[34] & B[354];
	B[173] = !B[173];
	B[355] = !B[3];
	B[174] = B[35] & B[355];
	B[356] = !B[174];
	B[175] = B[173] & B[356];
	B[357] = !B[33];
	B[176] = B[1] & B[357];
	B[176] = !B[176];
	B[358] = !B[32];
	B[177] = B[69] & B[358];
	B[178] = B[177] & B[0];
	B[178] = !B[178];
	B[179] = B[176] & B[178];
	B[179] = !B[179];
	B[180] = B[175] & B[179];
	B[180] = !B[180];
	B[359] = !B[35];
	B[181] = B[3] & B[359];
	B[181] = !B[181];
	B[182] = B[180] & B[181];
	B[183] = B[34] | B[174];
	B[183] = !B[183];
	B[184] = B[183] & B[2];
	B[184] = !B[184];
	B[185] = B[182] & B[184];
	B[185] = !B[185];
	B[360] = !B[167];
	B[186] = B[185] & B[360];
	B[187] = B[67] & B[161];
	B[188] = B[186] & B[187];
	B[188] = !B[188];
	B[189] = B[172] & B[188];
	B[190] = B[170] & B[189];
	B[190] = !B[190];
	B[191] = B[159] & B[190];
	B[192] = B[158] & B[191];
	B[192] = !B[192];
	B[193] = B[157] & B[192];
	B[194] = B[137] & B[193];
	B[194] = !B[194];
	B[361] = !B[113];
	B[195] = B[194] & B[361];
	B[362] = !B[120];
	B[196] = B[195] & B[362];
	B[363] = !B[123];
	B[197] = B[107] & B[363];
	B[198] = B[197] & B[68];
	B[199] = B[196] & B[198];
	B[199] = !B[199];
	B[200] = B[126] & B[199];
	B[201] = !B[200];
	B[202] = B[56] & B[55];
	B[203] = B[58] & B[57];
	B[204] = B[202] & B[203];
	B[205] = B[60] & B[59];
	B[206] = B[62] & B[61];
	B[207] = B[205] & B[206];
	B[208] = B[204] & B[207];
	B[208] = !B[208];
	B[364] = !B[82];
	B[209] = B[208] & B[364];
	B[209] = !B[209];
	B[365] = !B[209];
	B[210] = B[201] & B[365];
	B[211] = B[24] ^ B[56];
	B[212] = B[23] ^ B[85];
	B[212] = !B[212];
	B[213] = B[211] | B[212];
	B[213] = !B[213];
	B[214] = B[26] ^ B[58];
	B[215] = B[25] ^ B[57];
	B[216] = B[214] | B[215];
	B[216] = !B[216];
	B[217] = B[213] & B[216];
	B[218] = B[28] ^ B[60];
	B[219] = B[27] ^ B[59];
	B[220] = B[218] | B[219];
	B[220] = !B[220];
	B[221] = B[30] ^ B[62];
	B[222] = B[29] ^ B[61];
	B[223] = B[221] | B[222];
	B[223] = !B[223];
	B[224] = B[220] & B[223];
	B[225] = B[217] & B[224];
	B[225] = !B[225];
	B[366] = !B[225];
	B[226] = B[210] & B[366];
	B[226] = !B[226];
	B[227] = B[100] & B[226];
	B[228] = B[30] & B[94];
	B[228] = !B[228];
	B[229] = B[92] & B[27];
	B[229] = !B[229];
	B[230] = B[88] & B[25];
	B[230] = !B[230];
	B[231] = B[85] & B[23];
	B[231] = !B[231];
	B[367] = !B[231];
	B[232] = B[84] & B[367];
	B[232] = !B[232];
	B[368] = !B[24];
	B[233] = B[232] & B[368];
	B[233] = !B[233];
	B[234] = B[56] & B[231];
	B[234] = !B[234];
	B[235] = B[233] & B[234];
	B[236] = B[57] & B[72];
	B[236] = !B[236];
	B[237] = B[235] & B[236];
	B[237] = !B[237];
	B[238] = B[230] & B[237];
	B[239] = B[26] & B[87];
	B[239] = !B[239];
	B[240] = B[238] & B[239];
	B[240] = !B[240];
	B[241] = B[58] & B[71];
	B[241] = !B[241];
	B[242] = B[240] & B[241];
	B[243] = B[59] & B[76];
	B[243] = !B[243];
	B[244] = B[242] & B[243];
	B[244] = !B[244];
	B[245] = B[229] & B[244];
	B[246] = B[28] & B[91];
	B[246] = !B[246];
	B[247] = B[245] & B[246];
	B[247] = !B[247];
	B[248] = B[60] & B[75];
	B[248] = !B[248];
	B[249] = B[247] & B[248];
	B[250] = B[79] & B[61];
	B[250] = !B[250];
	B[251] = B[249] & B[250];
	B[251] = !B[251];
	B[252] = B[29] & B[95];
	B[252] = !B[252];
	B[253] = B[251] & B[252];
	B[253] = !B[253];
	B[254] = B[62] & B[78];
	B[254] = !B[254];
	B[255] = B[253] & B[254];
	B[255] = !B[255];
	B[256] = B[228] & B[255];
	B[256] = !B[256];
	B[369] = !B[256];
	B[257] = B[227] & B[369];
	B[257] = !B[257];
	B[258] = B[23] & B[24];
	B[259] = B[26] & B[25];
	B[260] = B[258] & B[259];
	B[261] = B[28] & B[27];
	B[262] = B[30] & B[29];
	B[263] = B[261] & B[262];
	B[264] = B[260] & B[263];
	B[264] = !B[264];
	B[265] = B[69] & B[173];
	B[370] = !B[0];
	B[266] = B[32] & B[370];
	B[266] = !B[266];
	B[267] = B[161] & B[266];
	B[268] = B[267] & B[67];
	B[269] = B[265] & B[268];
	B[371] = !B[3];
	B[270] = B[35] & B[371];
	B[270] = !B[270];
	B[271] = B[66] & B[270];
	B[272] = B[271] & B[166];
	B[273] = B[146] & B[165];
	B[274] = B[273] & B[145];
	B[275] = B[272] & B[274];
	B[276] = B[269] & B[275];
	B[277] = B[133] & B[141];
	B[278] = B[277] & B[129];
	B[279] = B[107] & B[64];
	B[280] = B[279] & B[68];
	B[281] = B[278] & B[280];
	B[282] = B[65] & B[134];
	B[283] = B[282] & B[112];
	B[372] = !B[22];
	B[284] = B[54] & B[372];
	B[284] = !B[284];
	B[285] = B[284] & B[111];
	B[286] = B[285] & B[101];
	B[287] = B[283] & B[286];
	B[288] = B[281] & B[287];
	B[289] = B[276] & B[288];
	B[289] = !B[289];
	B[290] = B[264] & B[289];
	B[373] = !B[201];
	B[291] = B[290] & B[373];
	B[292] = B[225] | B[291];
	B[293] = B[100] & B[292];
	B[374] = !B[256];
	B[294] = B[293] & B[374];
	B[295] = B[257] ^ B[294];
	B[296] = B[83] & B[295];
	B[297] = B[257] ^ B[296];
	B[375] = !B[291];
	B[298] = B[200] & B[375];
	B[298] = !B[298];
	B[376] = !B[209];
	B[299] = B[298] & B[376];
	B[299] = !B[299];
	B[300] = B[100] & B[299];
	B[377] = !B[225];
	B[301] = B[300] & B[377];
	B[302] = B[297] | B[301];
	B[302] = !B[302];
	Bit ret = B[302];

	delete[] B;

	return ret;	
}

