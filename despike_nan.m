function [y, locs] = despike_nan(x, smooth_win, thresh, diff_thresh, expand_win, varargin)
% DESPIKE_NAN
% スパイクノイズ区間を検出し，元波形の該当区間を NaN に置換する関数
%
% 構文:
%   y = despike_nan(x, smooth_win, thresh, diff_thresh, expand_win)
%   y = despike_nan(x, smooth_win, thresh, diff_thresh, expand_win, opts)
%   y = despike_nan(x, smooth_win, thresh, diff_thresh, expand_win, Name, Value, ...)
%   [y, locs] = despike_nan(___)
%
% 入力:
%   x           : 入力信号
%   smooth_win  : 低周波カット用の movmean 窓幅
%   thresh      : 低周波カット後の信号に対する閾値（NaN → 無効）
%   diff_thresh : 微分値に対する閾値（NaN → 無効）
%   expand_win  : マスク拡張用 movmean 窓幅
%                 スカラー: 対称窓  /  2要素ベクトル [nb nf]: 非対称窓
%                 NaN → 拡張なし
%
% オプション (Name-Value または構造体):
%   use_abs        : true のとき |x| に変換してから検出 (既定: false)
%   highpass_mode  : 'lowpass' | 'dc' | 'none'          (既定: 'lowpass')
%   lowpass_win    : 検出前の前平滑化窓幅（NaN = 無効）  (既定: NaN)
%   bilateral_mask : diff マスクを 1 サンプル前にも拡張  (既定: false)
%   min_spike_dur  : スパイクの最小連続サンプル数        (既定: 1)
%   max_spike_dur  : スパイクの最大連続サンプル数        (既定: Inf)
%   debug_flag     : true のとき検出過程をデバッグ表示   (既定: false)
%
% 出力:
%   y    : スパイク区間を NaN にした信号
%   locs : NaN を入れたサンプルのインデックス（列ベクトル）

    opts = parse_opts_(varargin);

    use_abs        = opts.use_abs;
    highpass_mode  = opts.highpass_mode;
    lowpass_win    = opts.lowpass_win;
    bilateral_mask = opts.bilateral_mask;
    min_spike_dur  = opts.min_spike_dur;
    max_spike_dur  = opts.max_spike_dur;
    debug_flag     = opts.debug_flag;

    was_row = isrow(x);
    x = x(:);

    if use_abs
        x_detect = abs(x);
    else
        x_detect = x;
    end

    [mask, x_hp_in, x_lf, x_hp, dx] = despike_mask_( ...
        x_detect, smooth_win, thresh, diff_thresh, expand_win, ...
        highpass_mode, lowpass_win, bilateral_mask);

    mask = duration_filter_(mask, min_spike_dur, max_spike_dur);

    y = x;
    y(mask) = NaN;

    if nargout >= 2
        locs = find(mask);
    end

    if debug_flag
        t = (1:length(x))';
        use_diff  = ~(isscalar(diff_thresh) && isnan(diff_thresh));
        n_panels  = 2 + use_diff;

        x_hp_in_valid = x_hp_in;
        x_hp_in_valid(mask) = NaN;

        figure;
        ax = gobjects(n_panels, 1);

        ax(1) = subplot(n_panels, 1, 1);
        plot(t, x_hp_in,       'b'); hold on;
        plot(t, x_lf,          'r', 'LineWidth', 1.2);
        plot(t, x_hp_in_valid, 'g', 'LineWidth', 1.2);
        legend('入力（ローパス後）', '低周波成分', '有効区間', 'Location', 'best');
        title(sprintf('低周波カット確認  [mode: %s, use_abs: %d]', highpass_mode, use_abs));

        ax(2) = subplot(n_panels, 1, 2);
        plot(t, x_hp, 'm'); hold on;
        if ~(isscalar(thresh) && isnan(thresh))
            yline( thresh, '--r');
            yline(-thresh, '--r');
        end
        title('低周波カット後（thresh 判定）');

        if use_diff
            ax(3) = subplot(n_panels, 1, 3);
            plot(t, dx, 'c'); hold on;
            yline( diff_thresh, '--r');
            yline(-diff_thresh, '--r');
            title('微分値（diff\_thresh 判定）');
        end

        linkaxes(ax, 'x');
    end

    if was_row
        y = y.';
    end

end


% =========================================================================
function opts = parse_opts_(args)
    opts.use_abs        = false;
    opts.highpass_mode  = 'lowpass';
    opts.lowpass_win    = NaN;
    opts.bilateral_mask = false;
    opts.min_spike_dur  = 1;
    opts.max_spike_dur  = Inf;
    opts.debug_flag     = false;

    if isempty(args)
        return;
    end

    idx = 1;

    if isstruct(args{1})
        s = args{1};
        fields = fieldnames(s);
        for k = 1:numel(fields)
            f = fields{k};
            if isfield(opts, f)
                opts.(f) = s.(f);
            else
                error('despike_nan: 未知のオプションフィールド ''%s''', f);
            end
        end
        idx = 2;
    elseif islogical(args{1}) || (isnumeric(args{1}) && isscalar(args{1}) && numel(args) == 1)
        opts.debug_flag = logical(args{1});
        return;
    end

    valid_names = fieldnames(opts);
    while idx <= numel(args)
        name = args{idx};
        if ~ischar(name) && ~isstring(name)
            error('despike_nan: オプションは ''Name'', Value の形式で指定してください（%d 番目の引数）', idx);
        end
        name = char(name);
        if ~ismember(name, valid_names)
            error('despike_nan: 未知のオプション名 ''%s''\n有効なオプション: %s', ...
                  name, strjoin(valid_names, ', '));
        end
        if idx + 1 > numel(args)
            error('despike_nan: オプション ''%s'' に対応する値がありません', name);
        end
        opts.(name) = args{idx + 1};
        idx = idx + 2;
    end
end


% =========================================================================
function [mask, x_hp_in, x_lf, x_hp, dx] = despike_mask_( ...
        x, smooth_win, thresh, diff_thresh, expand_win, ...
        highpass_mode, lowpass_win, bilateral_mask)

    if ~(isscalar(lowpass_win) && isnan(lowpass_win))
        x_hp_in = movmean(x, lowpass_win, 'omitnan');
    else
        x_hp_in = x;
    end

    switch lower(highpass_mode)
        case 'lowpass'
            x_lf = movmean(x_hp_in, smooth_win, 'omitnan');
            x_hp = x_hp_in - x_lf;
        case 'dc'
            x_lf = repmat(mean(x_hp_in, 'omitnan'), size(x_hp_in));
            x_hp = x_hp_in - x_lf;
        case 'none'
            x_lf = zeros(size(x_hp_in));
            x_hp = x_hp_in;
        otherwise
            error('despike_nan: highpass_mode は ''lowpass'', ''dc'', ''none'' のいずれかを指定してください。');
    end

    if isscalar(thresh) && isnan(thresh)
        mask_vel = false(size(x));
    else
        mask_vel = abs(x_hp) > thresh;
    end

    valid = ~isnan(x);
    if any(~valid) && sum(valid) >= 2
        t_all      = (1:length(x))';
        x_for_diff = interp1(t_all(valid), x(valid), t_all, 'spline');
    else
        x_for_diff = x;
    end
    dx = [NaN; diff(x_for_diff)];

    if isscalar(diff_thresh) && isnan(diff_thresh)
        mask_acc = false(size(x));
    else
        mask_acc = abs(dx) > diff_thresh;
        if bilateral_mask
            mask_acc = mask_acc | [false; mask_acc(1:end-1)];
        end
    end

    mask0 = mask_vel | mask_acc;

    if isscalar(expand_win) && isnan(expand_win)
        mask = mask0;
    else
        mask_expand = movmean(double(mask0), expand_win);
        mask = mask_expand ~= 0;
    end
end


% =========================================================================
function mask = duration_filter_(mask, min_dur, max_dur)
    if (min_dur <= 1) && isinf(max_dur)
        return;
    end

    d    = diff([0; mask(:); 0]);
    ons  = find(d ==  1);
    offs = find(d == -1) - 1;
    lens = offs - ons + 1;

    for k = 1:numel(ons)
        if lens(k) < min_dur || lens(k) > max_dur
            mask(ons(k):offs(k)) = false;
        end
    end
end
