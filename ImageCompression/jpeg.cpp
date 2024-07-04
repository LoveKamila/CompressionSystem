#include "PicReader.h"
#include "huffmancode.h"
#include<iostream>
#include <fstream>
#include <stdio.h>
#include<cmath>
#include<map>
#include<bitset>
using namespace std;

class jpegcompress {
    private:
        /*BYTE*&, 传入BYTE指针引用,用于接收像素信息
        * img_width 图像的宽度信息
        * img_height 图像的高度信息
          * 像素信息每四个一组(R G B A)    */
        PicReader imread;
        BYTE* data = nullptr;
        UINT img_width, img_height;
        
        //哈夫曼编码查找表
        map<unsigned char, string> mp_DC_lumin;  //亮度直流
        map<unsigned char, string> mp_AC_lumin;  //亮度交流
        map<unsigned char, string> mp_DC_chro;   //色度直流
        map<unsigned char, string> mp_AC_chro;   //色度交流

    public:
        string tmp = "";  //存放写入数据时的8位01串

        //读入图片
        void read_picture(char* argv[]);
        //文件头的写入
        void prewrite(ofstream &fout);
        //建立四张huffman编码查找表
        void huffmancode();
        //建立huffman编码查找表
        void set_huffmancode_list(int chooselist, const char* code_cnt, const unsigned char* code_value);
        //写入压缩后的图像信息
        void writecode(int zz[64], int chooselist, ofstream &fout);
        //得到该8*8小块的01串编码
        void getcode(int zz[64], int chooselist, string& result);
        //处理原图像素信息
        void process_pixel(ofstream &fout);
        //写入文件结束及释放空间
        void end(ofstream &fout);
};

int zz_Y[N * N] = { 0 }, zz_Cb[N * N] = { 0 }, zz_Cr[N * N] = { 0 };  //zigzag扫描排序后的数组
int preY = 0, preCb = 0, preCr = 0;  //存放上一小块第一个数据的原值（用于差分）

//读入图片
void jpegcompress::read_picture(char* argv[])
{
    imread.readPic(argv[2]);
    imread.getData(data, img_width, img_height);
}

//计算每小段RLE编码的第二个数对应的编码
string getcode_second(int val)
{
    string codeval = "";
    for (int k = abs(val); k > 0; k = k / 2) {
        codeval = codeval + (k % 2 == 1 ? '1' : '0');
    }

    if (val < 0)  //处理负数的情况
        for (int k = 0; k < (int)codeval.length(); k++)
            if (codeval[k] == '1')
                codeval[k] = '0';
            else
                codeval[k] = '1';

    //字符串倒序
    for (int i = 0; i < (int)codeval.length() / 2; i++) {
        char ch = codeval[i];
        codeval[i] = codeval[codeval.length() - i - 1];
        codeval[codeval.length() - i - 1] = ch;
    }
    return codeval;
}

static string to_binary_str(int code, int bitslen, string buf) {
    int mask = 1 << (bitslen - 1);
    for (int i = 0; i < bitslen; i++) {
        if (code & mask)
            buf += '1';
        else
            buf += '0';
        mask >>= 1;
    }
    buf[bitslen] = '0';
    return buf;
}

//建立huffman编码查找表
//chooselist区分直流交流亮度色度，code_cnt指向不同长度编码的个数数组，code_value为待编码的字符
void jpegcompress::set_huffmancode_list(int chooselist, const char* code_cnt, const unsigned char* code_value)
{
    int sym = 0, code = 0;
    string buf = "";
    for (int i = 0; i < 16; i++) {
        int num = (int)code_cnt[i];
        int bitslen = i + 1;
        for (int j = 0; j < num; j++) {
            unsigned char symbol = code_value[sym];
            string result = to_binary_str(code, bitslen, buf);

            if (chooselist == 1) 
                mp_DC_lumin.insert(make_pair(symbol, result));
            else if (chooselist == 2)
                mp_AC_lumin.insert(make_pair(symbol, result));
            else if (chooselist == 3)  
                mp_DC_chro.insert(make_pair(symbol, result));
            else  
                mp_AC_chro.insert(make_pair(symbol, result));
      
            code++;
            sym++;
        }
        code <<= 1;
    }
}

//建立四张huffman编码查找表
void jpegcompress::huffmancode()
{
    set_huffmancode_list(1, DC_lumin, DC_lumin_value);
    set_huffmancode_list(2, AC_lumin, AC_lumin_value);
    set_huffmancode_list(3, DC_chro, DC_chro_value);
    set_huffmancode_list(4, AC_chro, AC_chro_value);
}

//得到该8*8小块的01串编码
//zz为以zigzag方式排列成的一维数组，chooselist区分直流交流亮度色度，result储存编码
void jpegcompress::getcode(int zz[64], int chooselist, string& result)
{
    result = "";
    //numof0表示该小段最后一个数前面0的个数，val表示该段末尾的数
    int numof0 = 0, val = 0;  
    int end = 64;  //end表示EOB的位置
    for (int i = 63; i >= 0; i--)
        if (zz[i] != 0) {
            end = i + 1;
            break;
        }

    //DC部分
    string codeval = getcode_second(zz[0]);
    int len = codeval.length();
    if (len != 0) {
        if (chooselist == 1)
            result += mp_DC_lumin[(unsigned char)len];
        else
            result += mp_DC_chro[(unsigned char)len];
    }
    else {
        if (chooselist == 1)
            result += mp_DC_lumin[(unsigned char)0x00];
        else
            result += mp_DC_chro[(unsigned char)0x00];
    }
    result += codeval;

    //AC部分
    for (int i = 1; i < end; i++) {
        if (zz[i] == 0 && numof0 < 15)  //每小段以非0数结尾或者最多16个数
            numof0++;
        else {
            val = zz[i];
            string codeval = getcode_second(val);  //每小段RLE编码的第二个数对应的编码
            int len = codeval.length();  //BIT编码的中间一个数
            

            if (chooselist == 1) {
                    result += mp_AC_lumin[(unsigned char)((numof0 << 4) + len)];
            }
            else {
                    result += mp_AC_chro[(unsigned char)((numof0 << 4) + len)];
            }
            result += codeval;

            numof0 = 0;
        }
    }
}

//写入哈夫曼编码表
void print_huffman(char tag, const char  valuelen[], const unsigned char value[], ofstream &fout)
{
    fout.put(tag);  //表ID和表类型
    int sum = 0;    //数组value的元素个数为valuelen所有元素之和
    for (int i = 0; i < 16; i++) {
        fout.put(valuelen[i]);
        sum += (int)valuelen[i];
    }
    for (int i = 0; i < sum; i++)
        fout.put(value[i]);
}

//以蛇形输出量化表
void print_Q(const unsigned char Q[N][N], ofstream &fout)
{
    int tx = 0, ty = 0, dx = -1, dy = 1, cnt = 0;
    while (tx != N - 1 || ty != N - 1) {
        fout.put(Q[tx][ty]);
        if ((tx == 0 || tx == N - 1) && ty % 2 == 0) {
            ty++;
            dx *= -1;
            dy *= -1;
            continue;
        }
        if ((ty == 0 || ty == N - 1) && tx % 2 == 1) {
            tx++;
            dx *= -1;
            dy *= -1;
            continue;
        }
        tx += dx;
        ty += dy;
    }
    fout.put(Q[N - 1][N - 1]);
}

//文件头的写入
void jpegcompress::prewrite(ofstream &fout)
{
    /*SOI 图像开始*/
    fout.put((char)0xff), fout.put((char)0xd8);
    fout.put((char)0xff), fout.put((char)0xe0);  //APP0
    fout.put(0x00), fout.put(0x10);
    fout.put(0x4a), fout.put(0x46), fout.put(0x49), fout.put(0x46), fout.put(0x00);
    fout.put(0x01), fout.put(0x01), fout.put(0x00), fout.put(0x00), fout.put(0x01), fout.put(0x00), fout.put(0x01);
    fout.put(0x00), fout.put(0x00);

    /*DQT 量化表*/
    fout.put((char)0xff), fout.put((char)0xdb);  //标记代码
    fout.put((char)0x00), fout.put((char)0x84);  //数据长度(包括本字段)
    fout.put(0x00);  //精度及量化表ID
    print_Q(Qy, fout);    //以蛇形输出量化表Qy
    fout.put(0x01);  //精度及量化表ID
    print_Q(Qc, fout);    //以蛇形输出量化表Qc

    /*SOF0 帧图像开始*/
    fout.put((char)0xff), fout.put((char)0xc0);  //标记代码
    fout.put(0x00), fout.put(0x11), fout.put(0x08);
    fout.put(0x02); fout.put(0x00); //fout << img_height;  //图像高度(2字节)
    fout.put(0x02); fout.put(0x00); //fout << img_width;  //图像宽度(2字节)
    fout.put(0x03);  //颜色分量数
    //颜色分量信息
    fout.put(0x01), fout.put(0x11), fout.put(0x00); 
    fout.put(0x02), fout.put(0x11), fout.put(0x01);
    fout.put(0x03), fout.put(0x11), fout.put(0x01);

    /*DHT 定义哈夫曼表*/
    fout.put((char)0xff), fout.put((char)0xc4);  //标记代码
    fout.put((char)0x01), fout.put((char)0xa2);  //数据长度
    print_huffman(0x00, DC_lumin, DC_lumin_value, fout);
    print_huffman(0x10, AC_lumin, AC_lumin_value, fout);
    print_huffman(0x01, DC_chro, DC_chro_value, fout);
    print_huffman(0x11, AC_chro, AC_chro_value, fout);

    /*SOS 扫描开始*/
    fout.put((char)0xff), fout.put((char)0xda);  //标记代码
    fout.put(0x00), fout.put(0x0c);  //数据长度
    fout.put(0x03);  //颜色分量数
    fout.put(0x01), fout.put(0x00), fout.put(0x02), fout.put(0x11), fout.put(0x03), fout.put(0x11);  //颜色分量信息 (?)
    fout.put(0x00), fout.put(0x3f), fout.put(0x00);  //固定值
}

//写入压缩后的图像信息
//zz为以zigzag方式排列成的一维数组，chooselist区分直流交流亮度色度
void jpegcompress::writecode(int zz[64], int chooselist, ofstream &fout)
{
    string result = "";
    getcode(zz, chooselist, result);  //获取EOB前的编码，存入字符串result中
    if (chooselist == 1)
        result += mp_AC_lumin[(unsigned char)0x00];
    else
        result += mp_AC_chro[(unsigned char)0x00];
    //cout << result;
    //将result按单字节写入文件
    for (int pos = 0; pos < (int)result.length(); pos++) {
        tmp += result[pos];
        if ((int)tmp.length() >= 8) {   //tmp满8位则将该字节编码转为字符输出
            //cout << n << endl;
            bitset<8> n(tmp);
            unsigned char ch = (unsigned char)(n.to_ulong());
            fout.put(ch);
            if (ch == 0xff)
                fout.put(0x00);
            tmp = "";   //清空tmp
        }
    }
}

double alpha(int x)
{
    if (x == 0)
        return 1.0 / sqrt(8);
    else
        return 0.5;
}

//处理原图像素信息
void jpegcompress::process_pixel(ofstream &fout)
{
    int count = 0;  //当前第count个8*8小块
    for (DWORD i = 0; i < img_height; i += 8)
        for (DWORD j = 0; j < img_width; j += 8) {
            //图像分割：对每个8*8小块
            count++;
            double Y[N][N] = { 0 }, Cb[N][N] = { 0 }, Cr[N][N] = { 0 };
            for (DWORD u = 0; u < 8; u++)
                for (DWORD v = 0; v < 8; v++) {
                    int k = (i + u) * img_width * 4 + (j + v) * 4;
                    //颜色空间转换RGB->YCbCr
                    int R = data[k], G = data[k + 1], B = data[k + 2], A = data[k + 3];
                    Y[u][v] = 0.29871 * R + 0.58661 * G + 0.11448 * B - 128;
                    Cb[u][v] = -0.16874 * R - 0.33126 * G + 0.5 * B;
                    Cr[u][v] = 0.5 * R - 0.41869 * G - 0.08131 * B;
                }

            //离散余弦变换
            double f_Y[N][N] = { 0 }, f_Cb[N][N] = { 0 }, f_Cr[N][N] = { 0 };
            for (int u = 0; u < N; u++)
                for (int v = 0; v < N; v++) {
                    double tmp1 = 0, tmp2 = 0, tmp3 = 0;
                    for (int tu = 0; tu < N; tu++)
                        for (int tv = 0; tv < N; tv++) {
                            tmp1 += Y[tu][tv] * cos(1.0*(2 * tu + 1) * u * Pi / 16.0) * cos(1.0*(2 * tv + 1) * v * Pi / 16.0);
                            tmp2 += Cb[tu][tv] * cos(1.0*(2 * tu + 1) * u * Pi / 16.0) * cos(1.0*(2 * tv + 1) * v * Pi / 16.0);
                            tmp3 += Cr[tu][tv] * cos(1.0*(2 * tu + 1) * u * Pi / 16.0) * cos(1.0*(2 * tv + 1) * v * Pi / 16.0);
                        }
                    f_Y[u][v] = alpha(u) * alpha(v) * tmp1;
                    f_Cb[u][v] = alpha(u) * alpha(v) * tmp2;
                    f_Cr[u][v] = alpha(u) * alpha(v) * tmp3;
                }

            //数据量化
            int qua_Y[N][N], qua_Cb[N][N], qua_Cr[N][N];
            for (int u = 0; u < N; u++)
                for (int v = 0; v < N; v++) {
                    qua_Y[u][v] = (int)round(f_Y[u][v] / (int)Qy[u][v]);
                    qua_Cb[u][v] = (int)round(f_Cb[u][v] / (int)Qc[u][v]);
                    qua_Cr[u][v] = (int)round(f_Cr[u][v] / (int)Qc[u][v]);
                }
            
            //Zigzag扫描排序 (把量化后的二维矩阵转变成一个一维数组)

            int tmpY = qua_Y[0][0];
            int tmpCb = qua_Cb[0][0];
            int tmpCr = qua_Cr[0][0];
            zz_Y[0] = qua_Y[0][0] - preY;
            zz_Cb[0] = qua_Cb[0][0] - preCb;
            zz_Cr[0] = qua_Cr[0][0] - preCr;
            preY = tmpY;
            preCb = tmpCb;
            preCr = tmpCr;

            int tx = 0, ty = 1, dx = 1, dy = -1, cnt = 1;
            while (tx != N - 1 || ty != N - 1) {
                zz_Y[cnt] = qua_Y[tx][ty];
                zz_Cb[cnt] = qua_Cb[tx][ty];
                zz_Cr[cnt] = qua_Cr[tx][ty];
                cnt++;
                if ((tx == 0 || tx == N - 1) && ty % 2 == 0) {
                    ty++;
                    dx *= -1;
                    dy *= -1;
                    continue;
                }
                if ((ty == 0 || ty == N - 1) && tx % 2 == 1) {
                    tx++;
                    dx *= -1;
                    dy *= -1;
                    continue;
                }
                tx += dx;
                ty += dy;
            }
            zz_Y[cnt] = qua_Y[N - 1][N - 1];
            zz_Cb[cnt] = qua_Cb[N - 1][N - 1];
            zz_Cr[cnt] = qua_Cr[N - 1][N - 1];

            //写入压缩后的图像信息
            writecode(zz_Y, 1, fout);
            writecode(zz_Cb, 2, fout);
            writecode(zz_Cr, 2, fout);
        }
        if (tmp != "") {    //输出最后剩余不足8位的
            bitset<8> n(tmp);
            fout.put((unsigned char)n.to_ulong());
        }
}

//写入文件结束及释放空间
void jpegcompress::end(ofstream &fout)
{
    fout.put((char)0xff);
    fout.put((char)0xd9);
    delete[] data;
    data = nullptr;
}

int main(int argc, char* argv[])
{
    if (argc != 3) {
        cerr << "Please make sure the number of parameters is correct." << endl;
        return -1;
    }
    if (strcmp(argv[1], "-read") && strcmp(argv[1], "-compress")) {
        cerr << "Unknown parameter!\nCommand list:\nzip/unzip" << endl;
        return -1;
    }

    jpegcompress jpeg;

    if (strcmp(argv[1], "-compress") == 0) {
        //文件打开
        ifstream fin(argv[2], ios::binary); // 以二进制方式打开输入文件
        if (!fin) {  // 输出错误信息并退出
            cerr << "Can not open the input file!" << endl;
            return -1;
        }

        ofstream fout("lena.jpg", ios::binary);
        if (!fout) {
            cerr << "Can not open the output file!" << endl;
            return -1;
        }

        jpeg.read_picture(argv);
        jpeg.prewrite(fout);
        jpeg.huffmancode();
        jpeg.process_pixel(fout);
        jpeg.end(fout);
        cout << "complete!" << endl;

        fin.close();
        fout.close();
    }
    else if (strcmp(argv[1], "-read") == 0) {
        PicReader imread;
        BYTE* data = nullptr;
        UINT x, y;

        imread.readPic(argv[2]);
        imread.getData(data, x, y);
        imread.showPic(data, x, y);
        delete[] data;
    }
 
    return 0;
}