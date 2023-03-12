#pragma once

extern int FILTER_RADIUS;
extern float BLOOM_THRESHOLD;
extern float BLOOM_SCALE;
extern bool SHOW_ONLY_BLOOM_THRESHOLD;
extern bool SHOW_ONLY_BLOOM;
void applyBloomEffect(Screen& screen);
