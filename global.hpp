#pragma once

// グローバル変数宣言と初期化
//__declspec(selectany) double h = 29.0;
//__declspec(selectany) char g_iarr[] = { 1, 2, 3, 4 };
//__declspec(selectany) double g_darr[] = { 1.0, 2.0, 3.0, 4.0 };


// グローバル定数宣言と初期化
const double origX = 315.0;
const double origY = 253.7;
const double origZ = 130.0;
const double h = 29.0;
const char *filepath = "./../";
// const char *filepath = "./"; // for local test
const int PLpts = 1000;
const double R_target = 6.0;
const double R_MCP = 14./2.; // Y target
const double w_MCP = 419.; // Y target
const double R_FC = 10./2.; // SFC center
const double w_FC = 533.5; // SFC center
