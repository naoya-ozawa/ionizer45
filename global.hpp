#pragma once

// グローバル変数宣言と初期化
//__declspec(selectany) double h = 29.0;
//__declspec(selectany) char g_iarr[] = { 1, 2, 3, 4 };
//__declspec(selectany) double g_darr[] = { 1.0, 2.0, 3.0, 4.0 };


// グローバル定数宣言と初期化
const double h = 29.0;
const char *filepath = "./../";
// const char *filepath = "./"; // for local test
//const char *filepath = "./../SIMION_electrodes/online-setup/1.5mm_center_500V/online-1.5mmCentered-0.5kV-468V-";
const int PLpts = 1000;
const double R_target = 6.0;
const double R_MCP = 12.5;

// Set the origin point that was defined in the SIMION simulation parameters
const double SIMION_origin_x = 349.9;
const double SIMION_origin_y = 253.7;
const double SIMION_origin_z = 130.0;
