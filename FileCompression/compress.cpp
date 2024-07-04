#include <fstream>
#include <iostream>
#include<cstdio>
#include <cstdlib>
#include <string>
#include<vector>
#include<queue>
#include<map>
#include <bitset>
using namespace std;

//Huffman树的结点
struct node {
    long long freq;  //字符出现频数
    unsigned char chr;   //字符
    node* left;  //左孩子
    node* right;  //右孩子
    node* father;  //父节点
    node() {
        chr = 0;
        left = NULL;
        right = NULL;
        father = NULL;
        freq = -1;
    }
};

//优先队列维护结点，按各字符频数从小到大排序
struct cmp {
    bool operator()(const node& n1, const node& n2) {
        return n1.freq > n2.freq;
    }
};
priority_queue <node, vector<node>, cmp> que;

long long cnt[256] = { 0 };  //统计每个字符频数
string content;  //将待压缩的文件一次性读入字符串中

struct CODE {
    int len = 0;  //该字符的编码长度
    bitset<32> cd = 0;  //该字符的编码
};
struct CODE code[256];  //以字符的ascii码为下标，储存编码


//建立Huffman树
struct node* build_tree()
{
    //将出现过的字符插入优先队列
    for (int i = 0; i < 256; i++)
        if (cnt[i] > 0) {
            struct node* newnode = new(nothrow)(struct node);
            if (newnode == NULL) {
                cout << "wrong" << endl;
                exit(-1);
            }
            newnode->freq = cnt[i];
            newnode->chr = (unsigned char)i;
            que.push(*newnode);
        }
    //建立哈夫曼树
    while (!que.empty()) {
        struct node* top1 = new(nothrow)(struct node);
        struct node* top2 = new(nothrow)(struct node);
        struct node* fa = new(nothrow)(struct node);  //父节点
        if (top1 == NULL || top2 == NULL || fa == NULL) {
            cout << "wrong" << endl;
            exit(-1);
        }
        *top1 = que.top();
        que.pop();
        if (que.empty())
            return top1;
        *top2 = que.top();
        que.pop();
        top1->father = fa;
        top2->father = fa;
        fa->left = top1;
        fa->right = top2;  
        fa->freq = top1->freq + top2->freq;
        que.push(*fa);
    }
    return NULL;
}

//得到各字符的编码
void get_code
(struct node* root, int len, bitset<32> x)
{
    //到达叶子结点
    if (root->left == NULL && root->right == NULL) {
        int i = (int)root->chr;
        code[i].len = len;
        code[i].cd = x;
        return;
    }
    //向左递归
    x[len] = 0;
    get_code(root->left, len + 1, x);
    //向右递归
    x[len] = 1;
    get_code(root->right, len + 1, x);
    return;
}

//回收Huffman树空间
void cleartree(node* root)
{
    if (!root)
        return;
    cleartree(root->left);
    cleartree(root->right);
    delete(root);
}


int main(int argc, char* argv[]) {
    //带参主函数输入格式判断
    cout << "Zipper 0.001! Author: root" << endl;
    if (argc != 4) {
        cerr << "Please make sure the number of parameters is correct." << endl;
        return -1;
    }
    if (strcmp(argv[3], "zip") && strcmp(argv[3], "unzip")) {
        cerr << "Unknown parameter!\nCommand list:\nzip/unzip" << endl;
        return -1;
    }

    //文件打开
    ifstream fin(argv[1], ios::binary); // 以二进制方式打开输入文件
    if (!fin) {  // 输出错误信息并退出
        cerr << "Can not open the input file!" << endl;
        return -1;
    }
    ofstream fout(argv[2], ios::binary); // 以二进制方式打开输出文件
    if (!fout) {  // 输出错误信息并退出
        cerr << "Can not open the output file!" << endl;
        return -1;
    }

    //无损压缩与解压缩
    /*****************无损压缩过程******************/
    if (strcmp(argv[3], "zip") == 0) {
        istreambuf_iterator<char> beg(fin), end; 
        // 设置两个文件指针，指向开始和结束，以 char(一字节) 为步长
        string content(beg, end); 
        // 将文件全部读入 string 字符串

        //统计各字符频数
        int pos = 0;
        while (pos < (int)content.length()) {
            cnt[(int)content[pos]]++;
            pos++;
        }

        //建树与编码
        bitset<32> x = 0;  //储存每个字符的huffman编码
        struct node* root = build_tree();   //建立Huffman树
        get_code(root, 0, x);   //得到各字符编码
        
        //以8个二进制位为单位写入文件
        bitset<8> n = 0;  
        int bitnum = 0;  //当前二进制位存入bitset的位置
        pos = 0;  //原文件字符串下标
        while (pos < (int)content.length())
        {
            int asc = (int)content[pos];
            int i = code[asc].len;  //当前字符的01串编码长度
            for (int j = 0; j < i; j++) {
                n[bitnum] = code[asc].cd[j];
                bitnum++;
                if (bitnum >= 8) {  //写满8位就输出
                    fout.put((char)n.to_ulong());
                    n = 0;
                    bitnum = 0;
                }
            }
            pos++;
        }
        if (bitnum != 0)   //输出最后剩余不足8位的
            fout.put((char)n.to_ulong());

        //输出各字符频数用于解压缩时还原
        ofstream outtree("tree.log", ios::out);
        if (!outtree) {
            cerr << "Can not open the output file!" << endl;
            return -1;
        }
        for (int i = 0; i < 256; i++)
            outtree << cnt[i] << " ";
        outtree.close();

        //释放Huffman树空间
        cleartree(root);
    }

    /*****************解压缩过程******************/
    else if (strcmp(argv[3], "unzip") == 0) {
        //读入各字符的频数
        ifstream intree("tree.log", ios::in);
        if (!intree) {
            cerr << "Can not open the output file!" << endl;
            return -1;
        }
        for (int i = 0; i < 256; i++)
            intree >> cnt[i];
        intree.close();

        //重新建立Huffman树
        struct node* root = build_tree();
        if (root == NULL) {
            cerr << "wrong" << endl;
            exit(-1);
        }

        //读入压缩后的文件及还原操作
        istreambuf_iterator<char> beg(fin), end; // 设置两个文件指针，指向开始和结束，以 char(一字节) 为步长
        string s(beg, end); // 将文件全部读入 string 字符串
        int pos = 0;
        struct node* tmproot = root;  //表示当前节点
        while (pos < (int)s.length()) {
            bitset<8> x((unsigned long)s[pos]); //将当前字符转为01串
            int i = 0;
            while (i < 8) {
                //遇0则往左孩子，1则往右孩子
                if (x[i] == 0)
                    tmproot = tmproot->left;
                else
                    tmproot = tmproot->right;
                i++;
                //到达叶子结点就输出该路径所表示的字符
                if (tmproot->left == NULL && tmproot->right == NULL) {
                    fout.put(tmproot->chr);
                    tmproot = root;
                }
            }
            pos++;
        }

        //释放Huffman树空间
        cleartree(root);
    }

    // 全部操作完文件后关闭文件句柄
    fin.close();
    fout.close();

    cout << "Complete!" << endl;

    return 0;
}