#include "splender.h"
#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <iterator>
#include <limits>
#include <numeric>
#include <random>
#include <vector>

using namespace std;

const int deckn = 3, goal = 15;
const int FST = 0, TAKE = 1, BUY = 2, RES = 3;

namespace Gem {
const int all = 6, normal = 5, num = 4;
const int hold_max = 10, res_max = 3;
int remain[Gem::all];
}; // namespace Gem

class Board {
  public:
    static const int row = 4;
    card *operator[](int r) { return table[r]; }

  private:
    card table[deckn][row]; //,top[3][5];//0.9 0.7 0.5 0.3 0.1
};
Board table;

class Deck {
    vector<card> deck;
    int top_i;

  public:
    void init(vector &orig) {
        deck = move(orig);
        top_i = 4;
    }
    void pop() {
        if (top_i >= deck.size()) {
            cerr << "deck overpopped\n";
            exit(1);
        }
        ++top_i;
    }
    card operator[](int i) { // 0 for top, not vector's 0
        i += top_i;
        if (i >= deck.size() || i < top_i) {
            cerr << "invalid index\n";
            exit(1);
        }
        return deck[i];
    }
    bool empty() { return top_i == deck.size() ? true : false; }
};
Deck deck[3];

vector<card> deck[deckn];

struct profile {
    int gem[Gem::all];
    int bns[Gem::normal], bns_sum;
    int pts;
    int step[deckn][table.row];
    vector<card> res;
    void init() {
        fill(gem, gem + Gem::all, 0);
        fill(bns, bns + Gem::normal, 0);
        bns_sum = 0;
        pts = 0;
        res.clear();
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 4; ++j) step[i][j] = 0;
    }
};
profile my, op;
enum CardStat { INF = numeric_limits<int>::max() };

int stage;

// double demand[4]; //bonus to buy
//   when step not clear
//   base: lv1: step //lv2,3:card color!=0
//   mid: lv1,2:step //lv3:card color
//   last2: lv3:step lv1,lv2: card color

int cal_step(const card c, const profile &p) {
    int diff[4]{};
    for (int i = 0; i < 4; ++i) {
        if (c.cost[i] > Gem::num + p.bns[i]) {
            return INF;
        }
        else diff[i] = max(c.cost[i] - p.gem[i] - p.bns[i], 0);
    }
    return 0; // not complete

    // TODO compute least feasible steps (consider insufficiency)

    /*
    int diff[4]{},extra[4]{};
    bool unafford=0;
    for(int i=0;i<4;++i){
            unafford=c.gem[i]>gpos_max+p.bns[i];
            if(unafford){
                    extra[i]=c[i]-gpos_max-p.bns[i];
            }
            else diff[i]=max(c[i]-p.gem[i]-p.bns[i],0);
    }
    if(unafford&&c.score>=4&&){
            if(accumulate(extra,extra+4))
    }
    */
}

void init(vector<card> stack_1, vector<card> stack_2, vector<card> stack_3) {

    int cnt[deckn][Gem::normal]{};
    for (auto c : stack_1) {
        ++cnt[0][c.gem];
        cout << c.gem << ' ' << c.cost[0] << ' ' << c.cost[1] << ' ' << c.cost[2] << ' '
             << c.cost[3] << ' ' << c.cost[4] << '\n';
    }
    for (auto c : stack_2) {
        ++cnt[1][c.gem];
        cout << c.gem << ' ' << c.cost[0] << ' ' << c.cost[1] << ' ' << c.cost[2] << ' '
             << c.cost[3] << ' ' << c.cost[4] << '\n';
    }
    for (auto c : stack_3) {
        ++cnt[2][c.gem];
        cout << c.gem << ' ' << c.cost[0] << ' ' << c.cost[1] << ' ' << c.cost[2] << ' '
             << c.cost[3] << ' ' << c.cost[4] << '\n';
    }
    /*
for (int j = 0; j < deckn; ++j) {
    for (int i : cnt[j]) cout << i << ' ';
    cout << endl;
}
    */
    // init var
    my.init();
    op.init();
    stage = 1;
    // fill(demand, demand+5, 0);

    // set table and deck
    fill(Gem::remain, Gem::remain + Gem::all, 4);
    for (int i = 0; i < table.row; ++i) {
        table[0][i] = stack_1[i];
        table[1][i] = stack_2[i];
        table[2][i] = stack_3[i];
    }
    deck[0].init(stack_1); // stack moved
    deck[1].init(stack_2);
    deck[2].init(stack_3);

    // analyze table

    /*
for (int j = 0; j < table.row; ++j) {
    my.step[1][j] = cal_step(table[1][j]);
}
    */

    // analyze deck top
}

struct move player_move(struct move m) {
    // analyze table
    /*
for (int i = 0; i <= stage; ++i) {
    for (int j = 0; j < 4; ++j) {
        my.step[1][j]= cal_step(tab[1][j]);
    }
}
    */

    // formulate strategy

    // assign the strategy
    m.type = 1;
    m.gem[0] = 0;
    m.gem[1] = 0;
    m.gem[2] = 0;
    m.gem[3] = 0;
    m.gem[4] = 0;
    return m;
}
