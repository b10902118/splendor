#include <cassert>
#include <cstdio>
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include "splender.h"
#pragma GCC diagnostic pop
#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <iterator>
#include <limits>
#include <numeric>
#include <random>
#include <vector>

using namespace std;

constexpr int deckn = 3, goal = 15;

constexpr int depth[deckn]{15, 10, 7};
constexpr double rate[3]{0.96, 0.88, 0.7};

constexpr int FST = 0, TAKE = 1, BUY = 2, RES = 3;

namespace Gem {
constexpr int all = 6, normal = 5, num = 4;
constexpr int hold_max = 10, res_max = 3;
int remain[Gem::all];
}; // namespace Gem

class Board {
  public:
    static constexpr int row = 4;
    card *operator[](int r) { return table[r]; }

  private:
    card table[deckn][row]; //,top[3][5];//0.9 0.7 0.5 0.3 0.1
};
Board table;

class Deck {
    vector<card> deck;
    int top_i;

  public:
    void init(vector<card> &orig) {
        deck = std::move(orig);
        top_i = 4;
    }
    void pop() {
        if (top_i >= static_cast<int>(deck.size())) {
            cerr << "deck overpopped\n";
            exit(1);
        }
        ++top_i;
    }
    card operator[](int i) { // 0 for top, not vector's 0
        i += top_i;
        if (i >= static_cast<int>(deck.size()) || i < top_i) {
            cerr << "invalid index\n";
            exit(1);
        }
        return deck[i];
    }
    bool empty() {
        return top_i == static_cast<int>(deck.size()) ? true : false;
    }
};
Deck deck[3];

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

double supply[Gem::normal];
double demand[Gem::normal];
double value[deckn][table.row];

double fac[deckn][depth[0]];
void fac_init() {
    assert(depth[0] >= depth[1] && depth[0] >= depth[2]);
    for (int lv = 0; lv < deckn; ++lv) {
        for (int i = 0; i <= table.row; ++i)
            fac[lv][i] = 1; // first five
        for (int i = table.row + 1; i < depth[lv]; ++i)
            fac[lv][i] = fac[lv][i - 1] * rate[lv];
    }
}

// cal_step

double score1(const card &cd) {
    static constexpr double c[6] = {0,       1 * 1,    1.2 * 2,
                                    1.3 * 3, 1.45 * 4, 1.55 * 5};
    double score = 0;
    // cout << cd.gem << ": ";
    for (int i : cd.cost) {
        // cout << i << ' ';
        score += c[i];
    }
    // cout << 10 - score << endl;
    return 10 - score;
}

double score2(const card &cd) {
    static constexpr double c1 = 1.2, c2 = 1.3; // c1*n - c2*pt
    double score = 0;
    // cout << cd.gem << ": ";
    for (int i : cd.cost) {
        // cout << i << ' ';
        score += c1 * i;
    }
    // cout << score - c2 * cd.score << '\n';
    return 10 - (score - c2 * cd.score);
}

double score3(const card &cd) {
    static constexpr double c1 = 1.2, c2 = 1.3;
    double score = 0;
    for (int i : cd.cost) {
        score += c1 * (i - 1);
    }
    return 10 - (score - c2 * cd.score);
}

void init(vector<card> stack_1, vector<card> stack_2,
          vector<card> stack_3) {
    // init var

    my.init();
    op.init();
    stage = 1;
    fill(supply, supply + Gem::normal, 0);
    // fill(demand, demand+5, 0);

    // analyze here because index of deck start from 4
    // supply
    for (int i = 0; i < depth[0]; ++i) {
        supply[stack_1[i].gem] += score1(stack_1[i]);
        // cout << stack_1[i].gem << " += " << score1(stack_1[i])
        //<< '\n';
    }
    double sup1[Gem::normal] = {supply[0], supply[1], supply[2],
                                supply[3], supply[4]};
    for (int i = 0; i < depth[1]; ++i) {
        supply[stack_2[i].gem] += 0.5 * score2(stack_2[i]);
        // cout << stack_2[i].gem << " += " << score2(stack_2[i])
        //<< '\n';
    }

    for (int i = 0; i < Gem::normal; ++i)
        cout << "Gem-" << i << ": " << supply[i] << '\t' << sup1[i]
             << '\t' << supply[i] - sup1[i] << '\n';

    // demand
    for (int i = 0; i < depth[2]; ++i) {
    }

    // compute demand

    // set table and deck
    fill(Gem::remain, Gem::remain + Gem::all, 4);
    for (int i = 0; i < table.row; ++i) {
        table[0][i] = stack_1[i];
        table[1][i] = stack_2[i];
        table[2][i] = stack_3[i];
    }
    deck[0].init(stack_1); // stack "moved"
    deck[1].init(stack_2);
    deck[2].init(stack_3);

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
