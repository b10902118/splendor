#include "splender.h"
#include <algorithm>
#include <iostream>
#include <iterator>
#include <limits>
#include <numeric>
#include <vector>

using namespace std;

const int lvn = 3, goal = 15;
vector<card> deck[lvn] = {vector<card>(40), vector<card>(30),
                          vector<card>(20)};
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
    card table[lvn][row]; //,top[3][5];//0.9 0.7 0.5 0.3 0.1
};
Board table;

struct profile {
    int gem[Gem::all];
    int bns[Gem::normal], bns_sum;
    int pts;
    int step[lvn][table.row];
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

namespace Stage {
const int one = 1, two = 2, three = 3;
}
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

void init(vector<card> stack_1, vector<card> stack_2,
          vector<card> stack_3) {
    cout << stack_1.size() << ' ' << stack_2.size() << ' '
         << stack_3.size() << '\n';
    for (auto c : stack_1) {
        cout << c.score << ' ' << c.gem << ' ' << c.cost[0] << ' '
             << c.cost[1] << ' ' << c.cost[2] << ' ' << c.cost[3]
             << ' ' << c.cost[4] << '\n';
    }
    // init var
    my.init();
    op.init();
    stage = Stage::one;
    // fill(demand, demand+5, 0);

    // set table
    fill(Gem::remain, Gem::remain + Gem::all, 4);
    for (int i = 0; i < table.row; ++i)
        table[0][i] = std::move(stack_1[i]);
    for (int i = 0; i < table.row; ++i)
        table[1][i] = std::move(stack_2[i]);
    for (int i = 0; i < table.row; ++i)
        table[2][i] = std::move(stack_2[i]);
    // analyze table

    /*
for (int j = 0; j < table.row; ++j) {
    my.step[1][j] = cal_step(table[1][j]);
}
    */

    // analyze deck top

    // store deck
    for (int i = 0, n = deck[0].size(); i < n - table.row; ++i)
        deck[0][i] = std::move(stack_1[n - 1 - i]);
    for (int i = 0, n = deck[1].size(); i < n - table.row; ++i)
        deck[1][i] = std::move(stack_2[n - 1 - i]);
    for (int i = 0, n = deck[2].size(); i < n - table.row; ++i)
        deck[2][i] = std::move(stack_3[n - 1 - i]);
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
    m.gem[0] = 1;
    m.gem[1] = 1;
    m.gem[2] = 1;
    m.gem[3] = 0;
    m.gem[4] = 0;
    return m;
}
