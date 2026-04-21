clear; close all;

%% ============================================================
%% データ読み込み・セクション定義
%% ============================================================
load("Data1.mat");
fs = 1000;   % [Hz]

cal_section  = 1418110-fs : 1538110+fs;
data_section = 2718310-fs : 13590300+fs;   % OKS セット解析区間

oks     = -Data1_Ch4.values;
headcal = -Data1_Ch1.values(cal_section);
Reye    = -Data1_Ch3.values;
Leye    = -Data1_Ch6.values;

%% ============================================================
%% 前処理
%% ============================================================
Reye_s = movmean(Reye, 33);
Leye_s = movmean(Leye, 33);

ReyeV = my_bibun(Reye_s, fs);
LeyeV = my_bibun(Leye_s, fs);

%% ============================================================
%% キャリブレーション
%% ============================================================
ReyeV_cal_msk = despike_nan(ReyeV(cal_section), 201, 1.5, 1.2, 81);
LeyeV_cal_msk = despike_nan(LeyeV(cal_section), 201, 1.5, 1.2, 81);

headcal_cal = headcal * 36;   % [deg/s]

valid_r = ~isnan(headcal_cal) & ~isnan(ReyeV_cal_msk);
valid_l = ~isnan(headcal_cal) & ~isnan(LeyeV_cal_msk);

p_r = polyfit(headcal_cal(valid_r), ReyeV_cal_msk(valid_r), 1);
p_l = polyfit(headcal_cal(valid_l), LeyeV_cal_msk(valid_l), 1);

calib_gain = 1.17;
p_r(1) = p_r(1) * calib_gain;
p_l(1) = p_l(1) * calib_gain;

fprintf('p_r: slope=%.4f  intercept=%.4f\n', p_r(1), p_r(2));
fprintf('p_l: slope=%.4f  intercept=%.4f\n', p_l(1), p_l(2));

%% ============================================================
%% 全区間のサッカード・ノイズ除去 → deg/s 変換
%% ============================================================
ReyeV_msk = despike_nan(ReyeV, 201, 4, 0.5, [60 40], true);
LeyeV_msk = despike_nan(LeyeV, 201, 4, 0.5, [60 40]);

scale_r = abs(p_r(1));
scale_l = abs(p_l(1));
eye_sign_data = -1;

ReyeV_msk_deg_full = eye_sign_data * ReyeV_msk / scale_r;
LeyeV_msk_deg_full = eye_sign_data * LeyeV_msk / scale_l;

%% ============================================================
%% OKS 信号
%% ============================================================
oks_plot = oks * 20;
oks_Vt   = my_bibun(movmean(oks, 33), fs);

%% ============================================================
%% セット検出パラメータ  ← 必要に応じて変更
%% ============================================================
peaks_per_set   = 15;   % 1セット内の想定ピーク数
min_set_gap_s   = 150;  % セット間ギャップ閾値 [s]
min_peak_dist_s = 5;    % ピーク検出の最小間隔 [s]（セット内ピーク間）
M               = 20;   % ブロックサイズ（セット数）  ← 10 に変更可

%% ============================================================
%% エポック窓パラメータ（刺激開始 t = 0 を基準）
%%
%%   t = 0          : 各セットの最初の OKS ピーク（= 視覚刺激開始）
%%   t = -pre_s     : エポック開始（刺激開始の pre_s 秒前）
%%   t = +post_s    : エポック終了（刺激開始の post_s 秒後）
%%
%%   post_s は「次セットの刺激開始時刻 + post_extra 秒」として
%%   セット間隔から自動計算します．
%%   手動で固定値を使いたい場合は post_s_manual に秒数を入れてください．
%% ============================================================
pre_s        = 10;   % 刺激開始の何秒前からエポックを切るか [s]
post_extra   = 10;   % 次セットの刺激開始から何秒後まで含めるか [s]
post_s_manual = [];  % [] = 自動計算 / 数値を入れると固定値 [s]

%% ============================================================
%% OKS ピーク検出（方向自動判定）
%% ============================================================
oks_sec = oks_Vt(data_section);
vals    = oks_sec(~isnan(oks_sec));

if isempty(vals)
    error('data_section: OKS 速度信号がすべて NaN です。');
end

prom_th  = max(1.0, 2.0 * std(vals));
dist_smp = round(min_peak_dist_s * fs);
gap_smp  = round(min_set_gap_s   * fs);

[~, locs_pos] = findpeaks( oks_sec, 'MinPeakDistance', dist_smp, 'MinPeakProminence', prom_th);
[~, locs_neg] = findpeaks(-oks_sec, 'MinPeakDistance', dist_smp, 'MinPeakProminence', prom_th);

if numel(locs_pos) >= numel(locs_neg)
    all_locs  = locs_pos;
    direction = 'positive';
else
    all_locs  = locs_neg;
    direction = 'negative';
end

fprintf('OKS 方向: %s  (全ピーク数 = %d)\n', direction, numel(all_locs));

%% ============================================================
%% セット分割（連続ピーク間隔 > min_set_gap_s でセット境界）
%% ============================================================
if isempty(all_locs)
    error('ピークが検出されませんでした。min_peak_dist_s や prom_th を確認してください。');
end

is_new_set = [true; diff(all_locs(:)) > gap_smp];
set_id     = cumsum(is_new_set);
n_sets_raw = max(set_id);

onset_raw  = nan(n_sets_raw, 1);
offset_raw = nan(n_sets_raw, 1);
npeaks_raw = nan(n_sets_raw, 1);

for s = 1:n_sets_raw
    pk_in      = all_locs(set_id == s);
    onset_raw(s)  = pk_in(1);
    offset_raw(s) = pk_in(end);
    npeaks_raw(s) = numel(pk_in);
end

%% セット一覧をコンソール出力
fprintf('\n--- セット一覧 ---\n');
fprintf('%4s | %5s | %9s | %10s | %12s\n', 'Set', 'peaks', 'onset[s]', 'offset[s]', 'duration[s]');
for s = 1:n_sets_raw
    fprintf('%4d | %5d | %9.1f | %10.1f | %12.1f\n', ...
        s, npeaks_raw(s), onset_raw(s)/fs, offset_raw(s)/fs, ...
        (offset_raw(s) - onset_raw(s)) / fs);
end

%% ピーク数フィルタ（peaks_per_set の半数以上あれば有効）
valid_mask  = npeaks_raw >= round(peaks_per_set * 0.5);
onset_valid = onset_raw(valid_mask);
n_sets      = numel(onset_valid);

fprintf('\n有効セット数: %d / %d\n', n_sets, n_sets_raw);

if n_sets == 0
    error('有効なセットが見つかりません。peaks_per_set と min_set_gap_s を確認してください。');
end

%% ============================================================
%% セット間隔（onset-to-onset）から post_s を自動計算
%% ============================================================
if n_sets >= 2
    inter_onset_s = diff(onset_valid) / fs;
    iot_min    = min(inter_onset_s);
    iot_median = median(inter_onset_s);
    iot_max    = max(inter_onset_s);
    fprintf('\nセット開始間隔 (onset-to-onset):\n');
    fprintf('  min = %.1f s  /  median = %.1f s  /  max = %.1f s\n', ...
        iot_min, iot_median, iot_max);
else
    iot_min = [];
    fprintf('\nセット数が 1 のためセット間隔を計算できません。post_s_manual を指定してください。\n');
end

if ~isempty(post_s_manual)
    post_s = post_s_manual;
    fprintf('post_s: 手動設定 = %.1f s\n', post_s);
elseif ~isempty(iot_min)
    post_s = iot_min + post_extra;
    fprintf('post_s: 自動設定 = min(onset間隔) + post_extra = %.1f + %.1f = %.1f s\n', ...
        iot_min, post_extra, post_s);
else
    error('post_s を決定できません。post_s_manual に値を入れてください。');
end

fprintf('エポック窓: t = -%.1f s  ～  t = +%.1f s  (t=0 が刺激開始)\n', pre_s, post_s);

%% ============================================================
%% エポック切り出し
%% onset_valid は data_section 内のローカルインデックス
%% t0_data で全信号絶対インデックスへ変換
%% ============================================================
t0_data = data_section(1) - 1;

[epochs_eye, epochs_oks, t_epoch] = build_trigger_epochs( ...
    ReyeV_msk_deg_full, oks_plot, onset_valid, t0_data, fs, pre_s, post_s);

fprintf('エポック行列: %d サンプル × %d セット\n', size(epochs_eye, 1), size(epochs_eye, 2));

%% ============================================================
%% ブロック平均
%% ============================================================
stats = block_average_with_weights(epochs_eye, epochs_oks, M);

n_blocks      = stats.n_blocks;
add_ave       = stats.mean_sig;
add_ave_oks   = stats.mean_ref;
n_valid_block = stats.n_valid;
weight_block  = stats.weight;

fprintf('ブロック数: %d  (M = %d sets/block)\n', n_blocks, M);

if n_blocks == 0
    error('ブロック数が 0 です。n_sets=%d < M=%d を確認してください。', n_sets, M);
end

%% ============================================================
%% 全ブロック一覧表示
%% ============================================================
nCol     = min(4, n_blocks);
nRow     = ceil(n_blocks / nCol);
cmap_fig = turbo(256);
cmin_fig = 0;
cmax_fig = M;

fig1 = figure('Position', [50 50 min(1600, nCol*420) min(1000, nRow*340+120)]);
tl1  = tiledlayout(nRow, nCol, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl1, sprintf('OKS セット加算平均  (M = %d sets / block)', M), 'FontSize', 20);
xlabel(tl1, 'Time [s]', 'FontSize', 20);
ylabel(tl1, 'Eye Vel. [deg/s]', 'FontSize', 20);

for b = 1:n_blocks
    ax = nexttile; hold on;

    %% 個別セット散布
    Xb   = stats.blocks_sig{b};
    x_sc = repmat(t_epoch(:), 1, size(Xb, 2));
    vsc  = ~isnan(Xb);

    scatter(x_sc(vsc), Xb(vsc), 4, ...
        'MarkerEdgeColor', [0.78 0.78 0.78], ...
        'MarkerEdgeAlpha', 0.15, ...
        'HandleVisibility', 'off');

    %% OKS 参照信号
    plot(t_epoch, add_ave_oks(:,b), 'k', 'LineWidth', 0.8, 'DisplayName', 'OKS');

    %% 加算平均（有効セット数でカラーマップ）
    n_vld = n_valid_block(:,b);
    vline = ~isnan(t_epoch(:)) & ~isnan(add_ave(:,b));
    x_vl  = t_epoch(vline);
    y_vl  = add_ave(vline, b);
    c_vl  = n_vld(vline);

    surface([x_vl x_vl], [y_vl y_vl], zeros(sum(vline), 2), ...
        [c_vl c_vl], ...
        'FaceColor', 'none', 'EdgeColor', 'interp', 'LineWidth', 2.5, ...
        'DisplayName', '眼球速度');

    yline(0, '-', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.8, 'HandleVisibility', 'off');
    xline(0, ':', 'Color', [0.20 0.20 0.20], 'LineWidth', 1.0, 'HandleVisibility', 'off');
    if ~isempty(iot_min)
        xline(iot_min, '--', 'Color', [0.55 0.20 0.20], 'LineWidth', 0.8, ...
            'HandleVisibility', 'off');   % 次セット刺激開始
    end

    b_start = (b-1)*M + 1;
    b_end   = min(b*M, n_sets);
    title(sprintf('Block %d  (sets %d – %d)', b, b_start, b_end), 'FontSize', 14);

    colormap(ax, cmap_fig);
    caxis(ax, [cmin_fig cmax_fig]);
    set(ax, 'FontSize', 14, 'Box', 'off', 'TickDir', 'out', 'LineWidth', 0.8);
    hold off;
end

linkaxes(findobj(fig1, 'Type', 'axes'), 'xy');
xlim([-pre_s, post_s]);

%% colorbar（最初のタイルに付与）
cb = colorbar(nexttile(tl1, 1), 'eastoutside');
cb.Label.String = 'Valid sets';
cb.FontSize = 12;

%% ============================================================
%% ブロック重ね合わせ表示
%% ============================================================
colors_b = lines(n_blocks);

fig2 = figure('Position', [100 100 1000 600]);
hold on;

for b = 1:n_blocks
    b_start = (b-1)*M + 1;
    b_end   = min(b*M, n_sets);
    plot(t_epoch, add_ave(:,b), ...
        'Color', colors_b(b,:), ...
        'LineWidth', 1.5, ...
        'DisplayName', sprintf('Block %d  (sets %d–%d)', b, b_start, b_end));
end

plot(t_epoch, add_ave_oks(:,1), 'k--', 'LineWidth', 1.0, 'DisplayName', 'OKS');
yline(0, '-', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.8, 'HandleVisibility', 'off');
xline(0, ':',  'Color', [0.3 0.3 0.3], 'LineWidth', 1.0, 'HandleVisibility', 'off');
if ~isempty(iot_min)
    xline(iot_min, '--', 'Color', [0.55 0.20 0.20], 'LineWidth', 1.0, ...
        'Label', '次セット開始', 'HandleVisibility', 'off');
end

xlabel('Time [s]', 'FontSize', 18);
ylabel('Eye Vel. [deg/s]', 'FontSize', 18);
title(sprintf('ブロック別加算平均  (M = %d sets/block)  |  t=0: 刺激開始', M), 'FontSize', 18);
legend('Location', 'best', 'FontSize', 12);
set(gca, 'FontSize', 16, 'Box', 'off', 'TickDir', 'out', 'LineWidth', 0.8);
xlim([-pre_s, post_s]);
hold off;

%% ============================================================
%% 有効セット数の時間推移
%% ============================================================
fig3 = figure('Position', [150 150 900 400]);
hold on;

mean_n_valid = mean(n_valid_block, 1);   % ブロックごとの平均有効セット数
bar(1:n_blocks, mean_n_valid, 0.6, 'FaceColor', [0.20 0.45 0.75], 'EdgeColor', 'none');
yline(M, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8);

xlabel('Block', 'FontSize', 16);
ylabel('Mean valid sets', 'FontSize', 16);
title(sprintf('ブロックごとの平均有効セット数 (M = %d)', M), 'FontSize', 16);
xlim([0.4, n_blocks+0.6]);
ylim([0, M*1.1]);
set(gca, 'FontSize', 14, 'Box', 'off', 'TickDir', 'out', 'LineWidth', 0.8);
hold off;

%% ============================================================
%% local functions
%% ============================================================
function [epochs, ref_epochs, t_epoch] = build_trigger_epochs( ...
        signal, ref, locs, t0, fs, pre_s, post_s)
    pre_smp   = round(pre_s  * fs);
    post_smp  = round(post_s * fs);
    epoch_len = pre_smp + post_smp + 1;
    n_trials  = numel(locs);

    epochs     = NaN(epoch_len, n_trials);
    ref_epochs = NaN(epoch_len, n_trials);

    for ii = 1:n_trials
        g  = t0 + locs(ii);
        i1 = g - pre_smp;
        i2 = g + post_smp;

        if i1 < 1 || i2 > length(signal) || i2 > length(ref)
            continue;
        end

        epochs(:, ii)     = signal(i1:i2);
        ref_epochs(:, ii) = ref(i1:i2);
    end

    t_epoch = (-pre_smp : post_smp)' / fs;
end

function stats = block_average_with_weights(epochs, ref_epochs, M)
    [epoch_len, n_trials] = size(epochs);
    n_blocks = floor(n_trials / M);

    mean_sig = NaN(epoch_len, n_blocks);
    mean_ref = NaN(epoch_len, n_blocks);
    std_sig  = NaN(epoch_len, n_blocks);
    n_valid  = zeros(epoch_len, n_blocks);
    weight   = zeros(epoch_len, n_blocks);

    blocks_sig = cell(1, n_blocks);
    blocks_ref = cell(1, n_blocks);

    for b = 1:n_blocks
        idx = (b-1)*M + 1 : b*M;

        X  = epochs(:, idx);
        XR = ref_epochs(:, idx);

        blocks_sig{b} = X;
        blocks_ref{b} = XR;

        mean_sig(:,b) = mean(X,  2, 'omitnan');
        mean_ref(:,b) = mean(XR, 2, 'omitnan');
        std_sig(:,b)  = std(X,   0, 2, 'omitnan');
        n_valid(:,b)  = sum(~isnan(X), 2);
        weight(:,b)   = n_valid(:,b) / M;
    end

    stats.mean_sig   = mean_sig;
    stats.mean_ref   = mean_ref;
    stats.std_sig    = std_sig;
    stats.n_valid    = n_valid;
    stats.weight     = weight;
    stats.blocks_sig = blocks_sig;
    stats.blocks_ref = blocks_ref;
    stats.n_blocks   = n_blocks;
    stats.epoch_len  = epoch_len;
end
