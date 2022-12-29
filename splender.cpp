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

constexpr int deckn = 3, reserved = 3, goal = 15;

constexpr int depth[deckn]{15, 10, 7};
constexpr double rate[3]{0.96, 0.88, 0.7};
constexpr double buybase = 50, part = 0.4; // for the second deck

constexpr int TAKE = 1, BUY = 2, RES = 3;
constexpr int NOCARD = 0;

namespace Gem {
constexpr int all = 6, normal = 5, num = 4, joker = 5;
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
    card pop() {
        if (top_i >= static_cast<int>(deck.size())) {
            cerr << "deck empty\n";
            return (card){0, {0, 0, 0, 0, 0}, 0, 0};
        }
        return deck[top_i++];
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

struct ResCard {
    card c;
    int lv;
    ResCard(card c, int lv) : c(c), lv(lv) {}
};

struct profile {
    int gem[Gem::all], gem_sum;
    int bns[Gem::normal], bns_sum;
    int pts;
    int step[deckn][table.row];
    vector<ResCard> res;
    void init() {
        fill(gem, gem + Gem::all, 0);
        fill(bns, bns + Gem::normal, 0);
        bns_sum = 0;
        pts = 0;
        res.clear();
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 4; ++j) step[i][j] = 0;
    }
    void buy(card c) {
        for (int k = 0; k < Gem::normal; ++k) {
            gem[k] -= c.cost[k];
            Gem::remain[k] += min(c.cost[k], gem[k]);
            gem_sum -= c.cost[k];
            if (gem[k] < 0) {
                gem[Gem::joker] += gem[k];
                Gem::remain[Gem::joker] -= gem[k];
                gem[k] = 0;
            }
        }
        pts += c.score;
        ++bns[c.gem];
        ++bns_sum;
    }
    void resv(card c, int lv) {
        res.push_back(ResCard(c, lv));
        ++gem[Gem::joker];
        ++gem_sum;
    }
};
profile my, op;

int stage;

double demand[Gem::normal];

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
struct Step {
    int n;
    int pick[Gem::normal];
    int pickn;
};

Step cal_step(const card c, const profile &p) {
    int diff[Gem::normal]{}, ncolor = 0, match = 0, color_match = 0;
    Step ret = {0, {0, 0, 0, 0, 0}, 0};
    if (my.gem_sum <= 7) ret.pickn = 3;
    else if (my.gem_sum == 8) ret.pickn = 2;

    for (int i = 0; i < Gem::normal; ++i) {
        diff[i] = c.cost[i] - p.gem[i] - p.bns[i];
        if (diff[i] > 0) {
            ++ncolor;
            int nmatch = min(diff[i], Gem::remain[i]);
            if (ret.pickn == 2 && match < 2 && nmatch > match) {
                ret.pick[color_match] = 0;
                ret.pick[i] = min(2, nmatch);
                color_match = i;
            }
            else if (ret.pickn == 3 && match < 3 && nmatch > 0) {
                ++match;
                color_match = i;
                ++ret.pick[i];
            }
        }
        else if (diff[i] < 0) diff[i] = 0;
    }

    if (ret.pickn == 3 && match == 1) {
        if (diff[color_match] >= 2 && Gem::remain[color_match] >= 2) {
            ret.pickn = 2;
            ret.pick[color_match] = 2;
        }
    }
    // estimate best case: three and then two
    while (ncolor >= 2) {
        for (int i = 0, cnt = 0; i < Gem::normal && cnt < 3; ++i) {
            if (diff[i] > 0) {
                if (--diff[i] == 0) --ncolor;
                ++cnt;
            }
        }
        ++ret.n;
    }
    if (ncolor == 1) {
        for (int d : diff)
            if (d > 0) ret.n += (d & 1) ? d / 2 + 1 : d / 2;
    }
    return ret;
}

double score1(const card &cd) {
    static constexpr double c[6] = {0,       1 * 1,    1.2 * 2,
                                    1.3 * 3, 1.45 * 4, 1.55 * 5};
    double score = 0;
    // cout << cd.gem << ": ";
    for (int i = 0; i < Gem::normal; ++i) {
        // cout << cd.cost[i] << ' ';
        score += c[cd.cost[i] - my.bns[i]];
    }
    // cout << 10 - score << endl;
    return 10 - score;
}

double score2(const card &cd) {
    static constexpr double c1 = 1.2, c2 = 1.3; // c1*n - c2*pt
    double score = 0;
    // cout << cd.gem << ": ";
    for (int i = 0; i < Gem::normal; ++i) {
        int cost = cd.cost[i] - my.bns[i];
        // cout << cost<< ' ';
        if (cost > 0) score += c1 * cost;
    }
    // cout << score - c2 * cd.score << '\n';
    return 10 - (score - c2 * cd.score);
}

double score3(const card &cd) {
    static constexpr double c1 = 1.2, c2 = 1.3;
    double score = 0;
    cout << cd.gem << ": ";
    for (int i = 0; i < Gem::normal; ++i) {
        int cost = cd.cost[i] - my.bns[i];
        cout << cost << ' ';
        if (cost > 0) score += c1 * (cost - 1);
    }
    cout << score - c2 * cd.score << '\n';
    return 10 - (score - c2 * cd.score);
}

void init(vector<card> stack_1, vector<card> stack_2,
          vector<card> stack_3) {
    // init var
    fac_init();
    my.init();
    op.init();
    stage = 1;
    fill(demand, demand + 5, 0);

    // analyze here because index of deck start from 4
    // supply 1
    double supply[Gem::normal]{};
    for (int i = 0; i < depth[0]; ++i) {
        supply[stack_1[i].gem] += fac[0][i] * score1(stack_1[i]);
        // cout << stack_1[i].gem << " += " << score1(stack_1[i])
        //<< '\n';
    }
    double sup1[Gem::normal] = {supply[0], supply[1], supply[2],
                                supply[3], supply[4]};
    // supply and demand 2
    for (int i = 0; i < depth[1]; ++i) {
        double v = fac[1][i] * part * score2(stack_2[i]);
        supply[stack_2[i].gem] += v;
        for (int type = 0; type < Gem::normal; ++type) {
            demand[type] += v * stack_2[i].cost[type];
        }
        // cout << stack_2[i].gem << " += " << score2(stack_2[i])
        //<< '\n';
    }

    for (int i = 0; i < Gem::normal; ++i) {
        cout << "Gem-" << i << ": " << supply[i] << '\t' << sup1[i]
             << '\t' << supply[i] - sup1[i] << '\n';
    }

    cout << "demand before:\n";
    for (double i : demand) cout << i << ' ';
    cout << endl;
    // demand 3
    for (int i = 0; i < depth[2]; ++i) {
        double v = fac[2][i] * score3(stack_3[i]);
        for (int type = 0; type < Gem::normal; ++type) {
            cout << v * stack_3[i].cost[type] << ' ';
            demand[type] += v * stack_3[i].cost[type];
        }
        cout << endl << endl;
    }
    cout << "demand after:\n";
    for (double i : demand) cout << i << ' ';
    cout << endl;

    for (int i = 0; i < Gem::normal; ++i) demand[i] *= supply[i];
    cout << "demand times supply:\n";
    for (double i : demand) cout << i << ' ';
    cout << endl;
    // decide major and minor by supply and demand
    int fst = (demand[0] > demand[1]) ? 0 : 1, scd = 1 - fst;
    for (int i = 2; i < Gem::normal; ++i) {
        if (demand[i] > demand[fst]) {
            scd = fst;
            fst = i;
        }
        else if (demand[i] > demand[scd]) {
            scd = i;
        }
    }
    cout << fst << ' ' << scd << endl;
    for (int i = 0; i < Gem::normal; ++i) {
        if (i == fst) demand[i] = 10;
        else if (i == scd) demand[i] = 5;
        else demand[i] = 3;
        cout << demand[i] << ' ';
    }
    cout << endl;

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
}

pair<int, int> get_index(int cid) {
    for (int i = 0; i < deckn; ++i) {
        for (int j = 0; j < table.row; ++j) {
            if (table[i][j].id == cid) {
                return pair<int, int>(i, j);
            }
        }
    }
    for (int i = 0; i < my.res.size(); ++i) {
        if (my.res[i].c.id == cid) return pair<int, int>(reserved, i);
    }
    cerr << "cannot get index\n";
    exit(1);
}
bool iscard(card c) { return c.gem; }

void update(profile &p, pair<int, int> card_index, const int type) {
    if (card_index.first == reserved) { // buy reserved card
        assert(type == BUY);
        auto it = p.res.begin() + card_index.second;
        p.buy(it->c);
    }
    else {
        const int i = card_index.first, j = card_index.second;
        if (type == BUY) {
            p.buy(table[i][j]);
            table[i][j] = deck[i].pop();
        }
        else { // type==RES
            p.resv(table[i][j], i);
            table[i][j] = deck[i].pop();
        }
    }
}

bool affordable(profile &p, card c) {
    for (int i = 0, jk = 0; i < Gem::normal; ++i) {
        if (c.cost[i] > p.gem[i]) {
            jk += c.cost[i] - p.gem[i];
            if (jk > p.gem[Gem::joker]) return false;
        }
    }
    return true;
}

bool reservable(profile &p) {
    return p.res.size() < Gem::res_max && p.gem_sum < 10 &&
           Gem::remain[Gem::joker] > 0;
}

struct move buy(pair<int, int> card_index) {
    card c;
    if (card_index.first < reserved) {
        c = table[card_index.first][card_index.second];
    }
    else {
        c = my.res[card_index.second].c;
    }
    update(my, card_index, BUY);
    return (struct move){BUY, c.id, {0, 0, 0, 0, 0}};
}

struct move reserve(pair<int, int> card_index) {
    card c = table[card_index.first][card_index.second];
    update(my, card_index, RES);
    return (struct move){RES, c.id, {0, 0, 0, 0, 0}};
}

void print(profile &p) {
    cout << "gem: ";
    for (int i = 0; i < Gem::all; ++i) cout << p.gem[i] << ' ';
    cout << endl;
    cout << "bonus: ";
    for (int i = 0; i < Gem::normal; ++i) cout << p.bns[i] << ' ';
    cout << endl;
    cout << "points: " << p.pts << endl;
    cout << "card reserved:\n";
    for (auto cd : p.res) {
        cout << "Gem: " << cd.c.gem << " Cost: ";
        for (auto i : cd.c.cost) cout << i << ' ';
        cout << endl;
    }
    cout << endl;
}

struct move player_move(struct move m) {
    // record opponent's move
    if (m.type != 0) { // update board
        switch (m.type) {
        case TAKE:
            for (int i = 0; i < Gem::normal; ++i) {
                op.gem[i] += m.gem[i];
                Gem::remain[i] -= m.gem[i];
            }
            break;
        case BUY:
            update(op, get_index(m.card_id), BUY);
            break;
        case RES:
            update(op, get_index(m.card_id), RES);
            break;
        }
    }
    print(op);

    // compute value and step
    vector<pair<int, int>> card_index;
    // reserved cards
    double value[deckn + 1][table.row]{};
    for (int i = 0; i < my.res.size(); ++i) {
        if (my.res[i].lv == 0)
            value[reserved][i] = score1(my.res[i].c);
        else if (my.res[i].lv == 1)
            value[reserved][i] = score1(my.res[i].c);
        else value[reserved][i] = score1(my.res[i].c);
    }

    // table cards
    for (int i = 0; i < table.row; ++i) {
        if (iscard(table[0][i])) {
            value[0][i] = score1(table[0][i]);
            card_index.push_back(pair<int, int>(0, i));
        }
        if (iscard(table[1][i])) {
            value[1][i] = score2(table[1][i]);
            card_index.push_back(pair<int, int>(1, i));
        }
        if (iscard(table[2][i])) {
            value[2][i] = score3(table[2][i]);
            card_index.push_back(pair<int, int>(2, i));
        }
    }
    sort(card_index.begin(), card_index.end(),
         [value](pair<int, int> a, pair<int, int> b) -> bool {
             return value[a.first][a.second] >
                    value[b.first][b.second];
         }); // capture not sure

    // choose best move

    pair<int, int> bestres;
    int totake[5]{};
    bool totake_set = false;
    auto it = card_index.begin();
    for (; it != card_index.end() &&
           value[it->first][it->second] > buybase;
         ++it) {
        card c;
        if (it->first < reserved) {
            c = table[it->first][it->second];
        }
        else {
            c = my.res[it->second].c;
        }
        if (affordable(my, c)) return buy(*it);
        else {
            if (reservable(my)) {
                return reserve(*it);
            }
            else if (totake_set == false) {
                // TODO if not too far "cal_step"
                // totake=next_step() //need to compute available
                // pick
            }
        }
    }
    // if(can_pick)

    for (; it != card_index.end(); ++it) { // assume not reservable
        card c;
        if (it->first < reserved) {
            c = table[it->first][it->second];
        }
        else {
            c = my.res[it->second].c;
        }
        if (affordable(my, c)) return buy(*it);
    }

    // no valid operation
    m.type = 1;
    m.gem[0] = 1;
    m.gem[1] = 1;
    m.gem[2] = 1;
    m.gem[3] = 0;
    m.gem[4] = 0;
    return m;
}
