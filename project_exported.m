classdef project_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                      matlab.ui.Figure
        TabGroup                      matlab.ui.container.TabGroup
        NghiemTab                     matlab.ui.container.Tab
        SaisoEditField                matlab.ui.control.EditField
        SaisLabel                     matlab.ui.control.Label
        SolanEditField                matlab.ui.control.EditField
        SlnlpEditFieldLabel           matlab.ui.control.Label
        KetquaEditField               matlab.ui.control.EditField
        KtqunghimEditFieldLabel       matlab.ui.control.Label
        PhuongphapButtonGroup         matlab.ui.container.ButtonGroup
        DaycungButton                 matlab.ui.control.ToggleButton
        NewtonButton                  matlab.ui.control.ToggleButton
        LapButton                     matlab.ui.control.ToggleButton
        ChiadoiButton                 matlab.ui.control.ToggleButton
        PhanlyEditField               matlab.ui.control.EditField
        NhpkhongphnlynghimEditFieldLabel  matlab.ui.control.Label
        PhuongtrinhEditField          matlab.ui.control.EditField
        NhpphngtrnhEditFieldLabel     matlab.ui.control.Label
        UIAxes                        matlab.ui.control.UIAxes
        NoisuyTab                     matlab.ui.container.Tab
        KqdathucEditField             matlab.ui.control.EditField
        KtquathcnisuyEditFieldLabel   matlab.ui.control.Label
        KqnoisuyEditField             matlab.ui.control.NumericEditField
        KtqunisuyEditFieldLabel       matlab.ui.control.Label
        NhapnoisuyEditField           matlab.ui.control.NumericEditField
        NhpgitrcnnisuyEditFieldLabel  matlab.ui.control.Label
        NoisuyButtonGroup             matlab.ui.container.ButtonGroup
        LagrangeButton                matlab.ui.control.ToggleButton
        NewtonButton2                 matlab.ui.control.ToggleButton
        NhapyEditField                matlab.ui.control.EditField
        NhpdliuyEditFieldLabel        matlab.ui.control.Label
        NhapxEditField                matlab.ui.control.EditField
        NhpdliuxEditFieldLabel        matlab.ui.control.Label
        UIAxes2                       matlab.ui.control.UIAxes
        HoiquyTab                     matlab.ui.container.Tab
        KqdudoanEditField             matlab.ui.control.NumericEditField
        KtqudonEditFieldLabel         matlab.ui.control.Label
        NhapdudoanEditField           matlab.ui.control.NumericEditField
        NhpgitrcndonEditFieldLabel    matlab.ui.control.Label
        KqpthoiquyEditField           matlab.ui.control.EditField
        KtquphngtrnhhiquyEditFieldLabel  matlab.ui.control.Label
        HoiquyButtonGroup             matlab.ui.container.ButtonGroup
        LogaritButton                 matlab.ui.control.ToggleButton
        HammuButton                   matlab.ui.control.ToggleButton
        TuyentinhButton               matlab.ui.control.ToggleButton
        NhapyEditField2               matlab.ui.control.EditField
        NhpdliuyEditFieldLabel_2      matlab.ui.control.Label
        NhapxEditField2               matlab.ui.control.EditField
        NhpdliuxEditFieldLabel_2      matlab.ui.control.Label
        UIAxes3                       matlab.ui.control.UIAxes
        DaohamTab                     matlab.ui.container.Tab
        NhaphEditField                matlab.ui.control.NumericEditField
        NhpbchEditFieldLabel          matlab.ui.control.Label
        NhapdaohamEditField           matlab.ui.control.NumericEditField
        NhpgitrcntnhohmEditField_2Label  matlab.ui.control.Label
        KqdaohamEditField             matlab.ui.control.NumericEditField
        KtquEditFieldLabel            matlab.ui.control.Label
        ChondaohamButtonGroup         matlab.ui.container.ButtonGroup
        TrungtmButton                 matlab.ui.control.ToggleButton
        LiButton                      matlab.ui.control.ToggleButton
        XpxtinButton                  matlab.ui.control.ToggleButton
        ChonsaisoButtonGroup          matlab.ui.container.ButtonGroup
        Oh2Button                     matlab.ui.control.ToggleButton
        OhButton                      matlab.ui.control.ToggleButton
        NhaphamEditField              matlab.ui.control.EditField
        NhphmsEditFieldLabel          matlab.ui.control.Label
        HocLabel                      matlab.ui.control.Label
        NhapyEditField3               matlab.ui.control.EditField
        NhpdliuyEditFieldLabel_3      matlab.ui.control.Label
        NhapxEditField3               matlab.ui.control.EditField
        NhpdliuxEditFieldLabel_3      matlab.ui.control.Label
        TichphanTab                   matlab.ui.container.Tab
        NhapcanEditField              matlab.ui.control.EditField
        NhpcntnhtchphnEditFieldLabel  matlab.ui.control.Label
        KqtichphanEditField           matlab.ui.control.NumericEditField
        KtquEditFieldLabel_2          matlab.ui.control.Label
        NhapNEditField                matlab.ui.control.NumericEditField
        NhpNEditFieldLabel            matlab.ui.control.Label
        ChontichphanButtonGroup       matlab.ui.container.ButtonGroup
        Simpson38Button               matlab.ui.control.ToggleButton
        Simpson13Button               matlab.ui.control.ToggleButton
        HnhthangButton                matlab.ui.control.ToggleButton
        Nhapham2EditField             matlab.ui.control.EditField
        NhphmsEditFieldLabel_2        matlab.ui.control.Label
        HocLabel_2                    matlab.ui.control.Label
        NhapyEditField4               matlab.ui.control.EditField
        NhpdliuyEditFieldLabel_4      matlab.ui.control.Label
        NhapxEditField4               matlab.ui.control.EditField
        NhpdliuxEditFieldLabel_4      matlab.ui.control.Label
        GioithieunhomTab              matlab.ui.container.Tab
        nthitkapptrnMatlabvccphntnhtonhcnhNghimNisuyHiquyohmTchphnLabel  matlab.ui.control.Label
        GiithiuLabel                  matlab.ui.control.Label
        ThitkvcodetonbnthchnhphngphptnhLabel  matlab.ui.control.Label
        NhimvLabel                    matlab.ui.control.Label
        NguynQucTrng22200172Label     matlab.ui.control.Label
        ThnhvinLabel                  matlab.ui.control.Label
        Nhm7Label                     matlab.ui.control.Label
    end

methods (Access = private)
    %Nghiem
    % Phương pháp chia đôi
    function [c, n] = chiadoi(app, fx, a, b, saiso)
        n = 0; % Số lần lặp
        e = abs(b - a);
        while e >= saiso
            c = (a + b) / 2;
            if fx(c) * fx(a) < 0
                b = c;
            else
                a = c;
            end
            e = abs(b - a);
            n = n + 1;
        end
    end
    % Phương pháp lặp
    function [c, n] = lap(app, fx, a, b, saiso)
        syms x;
        df = matlabFunction(diff(fx(x))); % Đạo hàm fx(x)
        fp = @(x) x - fx(x) / df(x); % Sinh hàm fp(x)
        x0 = (a + b) / 2; % Giá trị khởi tạo
        x1 = fp(x0);
        n = 0; % Số lần lặp
        while abs(x1 - x0) >= saiso
            x0 = x1;
            n = n + 1;
            x1 = fp(x0);
        end
        c = x1;
    end
    % Phương pháp tiếp tuyến (Newton)
    function [c, n] = newton(app, fx, a, b, saiso)
        syms x;
        df = matlabFunction(diff(fx(x))); % Đạo hàm fx(x)
        x0 = (a + b) / 2; % Giá trị khởi tạo
        x1 = x0 - fx(x0) / df(x0);
        n = 0; % Số lần lặp
        while abs(x1 - x0) >= saiso
            x0 = x1;
            n = n + 1;
            x1 = x0 - fx(x0) / df(x0);
        end
        c = x1;
    end
    % Phương pháp dây cung
    function [c, n] = daycung(app, fx, a, b, saiso)
        x0 = a;
        x1 = b;
        e = abs(x1 - x0);
        n = 0; % Số lần lặp
        while e >= saiso
            x0 = x1;
            x1 = (a * fx(b) - b * fx(a)) / (fx(b) - fx(a)); % Công thức dây cung
            if fx(x1) * fx(a) < 0
                b = x1;
            else
                a = x1;
            end
            n = n + 1;
            e = abs(x1 - x0);
        end
        c = x1;
    end
    %Noi suy
    %Đa thức Lagrange
    function result = Lagrange(app,xa,ya,x)
        n = length(xa);
        sum = 0;
        for i = 1:n
            product = ya(i);
            for j = 1:n
                if i~=j
                    product = product*(x-xa(j))/(xa(i) - xa(j));
                end
            end
            sum = sum + product;
        end
        result = sum;
    end
    function L = LagrangePolynomial(app,xa,ya)
        syms x;
        n = length(xa);
        L = 0;
        for i = 1:n
            Li = 1;
            for j = 1:n
                if i ~= j
                    Li = Li * (x - xa(j)) / (xa(i) - xa(j));
                end
            end
            L = L + Li * ya(i);
        end
        L = expand(L);
    end
    %Đa thức Newton
    function da = DividedDifference(app,xa, ya)
        n = length(xa);
        da = ya;
        for i = 1:n
            for j = 1:i-1
            da(i) = (da(j)-da(i))/(xa(j)-xa(i));
            end
        end
    end
    function result = NewtonForm(app,xa, da, x)
        n = length(da);
        result = da(n);
        for i = n-1:-1:1
            result = result*(x-xa(i))+da(i);
        end
    end
    function result = NewtonInterpolation(app,xa, ya, x)
        da = app.DividedDifference(xa,ya);
        result = app.NewtonForm(xa,da,x);
    end
    function polynomial = NewtonPolynomial(app,xa, ya)
        syms x;
        da = app.DividedDifference(xa, ya);
        n = length(xa);
        polynomial = da(1);
        term = 1;
        for i = 2:n
            term = term * (x - xa(i-1));
            polynomial = polynomial + da(i) * term;
        end
        polynomial = expand(polynomial);
    end
    function [a1, a0, r2] = Hoiquytuyentinh(app,x, y)
        n = length(x);
        sumx = 0;
        sumy = 0;       
        sumxy = 0;
        sumx2 = 0;
        st = 0;
        sr = 0;
        for i = 1:n
            sumx = sumx + x(i);
            sumy = sumy + y(i);
            sumxy = sumxy + x(i)*y(i);
            sumx2 = sumx2 + x(i)*x(i);
        end
        xm = sumx/n;
        ym = sumy/n;
        a1 = (n*sumxy - sumx*sumy)/(n*sumx2 - sumx*sumx);
        a0 = ym - a1*xm;
        for i = 1:n
            st = st + (y(i)- ym)^2;
            sr = sr + (y(i) - a1*x(i) - a0)^2;
        end
        r2 = (st - sr)/st;
    end
    function [a, b, r2] = Hoiquyhammu(app, x, y)
    % Hồi quy hàm mũ y = a * exp(b * x)
    % Chuyển về tuyến tính ln(y) = ln(a) + b * x
    if any(y <= 0)
        error('Dữ liệu y phải lớn hơn 0 để thực hiện hồi quy hàm mũ.');
    end    
    Y = log(y);  % Logarit hóa y
    [b, ln_a, r2] = Hoiquytuyentinh(app, x, Y); % Hồi quy tuyến tính trên x và log(y)
    a = exp(ln_a);  % Lấy log(a) thành a
    end
    function [a, b, r2] = Hoiquyhamlogarit(app, x, y)
    % Hồi quy hàm logarit y = a * ln(x) + b
    % x phải lớn hơn 0
    if any(x <= 0)
        error('Dữ liệu x phải lớn hơn 0 để thực hiện hồi quy hàm logarit.');
    end
    X = log(x);  % Logarit hóa x
    [a, b, r2] = Hoiquytuyentinh(app, X, y);  % Hồi quy tuyến tính trên log(x) và y
    end
    function I = tichphanhinhthang(app,fx,a,b,N)
        h = (b - a)/N;
        x_values = a:h:b;
        y_values = fx(x_values);
        I = h*(2*sum(y_values) - (y_values(1) + y_values(end)))/2;
    end
    function I = tichphanSimpson1_3(app,fx,a,b,N)
        if mod(N, 2) ~= 0
            error('Số đoạn chia N phải là số chẵn!')
        end
        h = (b - a)/N;
        x_values = a:h:b;
        y_values = fx(x_values);
        I = (h/3)*(y_values(1) + y_values(end) + 4*sum(y_values(2:2:end-1)) + 2*sum(y_values(3:2:end-2)));
    end
    function I = tichphanSimpson3_8(app,fx,a,b,N)
        if mod(N,3) ~= 0
            error('Số đoạn chia N phải là bội số của 3!');
        end
        h = (b - a)/N;
        x_values = a:h:b;
        y_values = fx(x_values);
        I = (3/8)*h*(y_values(1) + y_values(end) + 3*sum(y_values(2:3:end)) + 3*sum(y_values(3:3:end-1)) + 2*sum(y_values(4:3:end-3)));
    end    
end   

    % Callbacks that handle component events
    methods (Access = private)

        % Selection changed function: PhuongphapButtonGroup
        function PhuongphapButtonGroupSelectionChanged(app, event)
            selectedButton = app.PhuongphapButtonGroup.SelectedObject;
            try
                % Xử lý đầu vào
                fx = str2func(['@(x)', app.PhuongtrinhEditField.Value]); % Hàm số
                phanly = str2num(app.PhanlyEditField.Value); % Chuyển khoảng phân ly
                a = phanly(1);
                b = phanly(2);
                saiso = str2double(app.SaisoEditField.Value); % Chuyển sai số
                % Tìm nghiệm
                c = 0;
                n = 0;
                switch selectedButton.Text
                   case 'Chia đôi'
                        [c, n] = app.chiadoi(fx, a, b, saiso);
                   case 'Lặp'
                        [c, n] = app.lap(fx, a, b, saiso);
                   case 'Newton (tiếp tuyến)'
                        [c, n] = app.newton(fx, a, b, saiso);
                   case 'Dây cung'
                        [c, n] = app.daycung(fx, a, b, saiso);
                end
                % Cập nhật kết quả
                app.KetquaEditField.Value = num2str(c);
                app.SolanEditField.Value = num2str(n);
                % Vẽ đồ thị hàm số
                x = linspace(a, b, 500); % Tạo 500 điểm trong khoảng [a, b]
                y = arrayfun(fx, x); % Tính giá trị y = f(x) tại các điểm x
                plot(app.UIAxes, x, y, 'b-', 'LineWidth', 1.5); % Vẽ đồ thị
                hold(app.UIAxes, 'on'); % Cho phép vẽ thêm
                plot(app.UIAxes, c, fx(c), 'rx', 'MarkerSize', 10, 'LineWidth', 2); % Vẽ nghiệm
                hold(app.UIAxes, 'off'); % Tắt chế độ vẽ thêm
            catch ME
                % Hiển thị thông báo lỗi
                uialert(app.UIFigure, 'Có lỗi xảy ra! Kiểm tra lại đầu vào.', 'Lỗi');
                disp(ME.message); % Hiển thị lỗi trong Command Window
            end
        end

        % Selection changed function: NoisuyButtonGroup
        function NoisuyButtonGroupSelectionChanged(app, event)
            selectedButton = app.NoisuyButtonGroup.SelectedObject;
            try
                xa = str2num(app.NhapxEditField.Value);
                ya = str2num(app.NhapyEditField.Value);
                x = app.NhapnoisuyEditField.Value;
                if isempty(xa)||isempty(ya)||isempty(x)
                    error('Dữ liệu nhập không hợp lệ hoặc bị thiếu!');
                end
                if length(xa)~=length(ya)
                    error('Dữ liệu x và y phải có cùng số lượng phần tử!');
                end
                result = 0;
                polynomial = '';
                switch selectedButton.Text
                    case 'Lagrange'
                        result = app.Lagrange(xa,ya,x);
                        polynomial = app.LagrangePolynomial(xa,ya);
                    case 'Newton'
                        result = app.NewtonInterpolation(xa,ya,x);
                        polynomial = app.NewtonPolynomial(xa,ya);
                end
                app.KqnoisuyEditField.Value = result;
                app.KqdathucEditField.Value = char(polynomial);
                xi = linspace(min(xa), max(xa), 500); % Tạo các điểm trong khoảng [min(xa), max(xa)]
                switch selectedButton.Text
                    case 'Lagrange'
                        yi = arrayfun(@(t) app.Lagrange(xa, ya, t), xi); % Dùng Lagrange
                    case 'Newton'
                        yi = arrayfun(@(t) app.NewtonInterpolation(xa, ya, t), xi); % Dùng Newton
                end
                plot(app.UIAxes2, xa, ya, 'ro', 'MarkerSize', 8, 'LineWidth', 1.5); % Vẽ dữ liệu thực
                hold(app.UIAxes2, 'on');
                plot(app.UIAxes2, xi, yi, 'b-', 'LineWidth', 1.5); % Vẽ đồ thị nội suy
                plot(app.UIAxes2, x, result, 'gx', 'MarkerSize', 10, 'LineWidth', 2); % Vẽ điểm nội suy
                hold(app.UIAxes2, 'off');
                legend(app.UIAxes2, 'Dữ liệu thực', 'Đường nội suy', 'Điểm nội suy');
            catch ME
                % Hiển thị thông báo lỗi nếu có
                uialert(app.UIFigure, 'Có lỗi xảy ra! Vui lòng kiểm tra lại dữ liệu.', 'Lỗi');
                disp(ME.message); % Ghi lỗi vào Command Window để dễ dàng kiểm tra
            end
        end

        % Selection changed function: HoiquyButtonGroup
        function HoiquyButtonGroupSelectionChanged(app, event)
            selectedButton = app.HoiquyButtonGroup.SelectedObject;   
            try
                xi = str2num(app.NhapxEditField2.Value);
                yi = str2num(app.NhapyEditField2.Value);
                xdudoan = app.NhapdudoanEditField.Value;
                if isempty(xi) || isempty(yi) || isempty(xdudoan)
                    error('Dữ liệu nhập không hợp lệ hoặc bị thiếu!');
                end
                if length(xi)~=length(yi)
                    error('Dữ liệu x và y phải có cùng số lượng phần tử!');
                end
                pthoiquy = '';
                r2 = 0;
                yhoiquy = []; %hàm hồi quy
                switch selectedButton.Text
                    case 'Tuyến tính'
                        [a,b,r2] = app.Hoiquytuyentinh(xi,yi);
                        pthoiquy = sprintf('y = %.4fx + %.4f',a,b);
                        yhoiquy = @(x) a*x + b;
                    case 'Hàm mũ'
                        [a,b,r2] = app.Hoiquyhammu(xi,yi);
                        pthoiquy = sprintf('y = %.4fe^{%.4fx}',a,b);
                        yhoiquy = @(x) a*exp(b*x);
                    case 'Logarit'
                        [a,b,r2] = app.Hoiquyhamlogarit(xi,yi);
                        pthoiquy = sprintf('y = %.4fln(x) + %.4f',a,b);
                        yhoiquy = @(x) a*log(x) + b;
                end
                app.KqpthoiquyEditField.Value = pthoiquy;
                ydudoan = yhoiquy(xdudoan);
                app.KqdudoanEditField.Value = ydudoan;
                cla(app.UIAxes3);
                plot(app.UIAxes3,xi,yi,'bo', 'DisplayName','Dữ liệu thực');
                hold (app.UIAxes3, 'on');
                xplot = linspace(min(xi),max(xi),100); %tạo dải giá trị x
                plot(app.UIAxes3,xplot,yhoiquy(xplot), 'r-', 'DisplayName', 'Hàm hồi quy');
                legend(app.UIAxes3, 'show');
                hold (app.UIAxes3, 'off');  
            catch ME
                % Hiển thị thông báo lỗi nếu có
                uialert(app.UIFigure, 'Có lỗi xảy ra! Vui lòng kiểm tra lại dữ liệu.', 'Lỗi');
                disp(ME.message); % Ghi lỗi vào Command Window để dễ dàng kiểm tra
            end
        end

        % Selection changed function: ChondaohamButtonGroup, 
        % ...and 1 other component
        function ChondaohamButtonGroupSelectionChanged(app, event)
            selectedButton = app.ChondaohamButtonGroup.SelectedObject;
            try
                selectedButton2 = app.ChonsaisoButtonGroup.SelectedObject;
                xi = str2num(app.NhapxEditField3.Value);
                yi = str2num(app.NhapyEditField3.Value);
                func = app.NhaphamEditField.Value;
                h = app.NhaphEditField.Value;
                xdaoham = app.NhapdaohamEditField.Value;
                if isempty(xdaoham) || isnan(xdaoham)
                    error('Hãy nhập giá trị cần tính đạo hàm!');
                end
                if ~isempty(xi) && ~isempty(yi)
                    if length(xi) ~= length(yi)
                        error('Dữ liệu x và y phải có cùng số lượng phần tử!');
                    end
                    % Tìm hai điểm gần nhất với x_tinh để nội suy đạo hàm
                    [~, idx] = sort(abs(xi - xdaoham));
                    if length(idx) < 2
                        error('Dữ liệu không đủ để tính đạo hàm!');
                    end
                    x1 = xi(idx(1));
                    x2 = xi(idx(2));
                    y1 = yi(idx(1));
                    y2 = yi(idx(2));
                    % Tính đạo hàm theo phương pháp
                    switch selectedButton.Text
                        case 'Xấp xỉ tiến'
                            result = (y2 - y1) / (x2 - x1);
                        case 'Lùi'
                            result = (y2 - y1) / (x2 - x1);
                        case 'Trung tâm'
                            if length(xi) < 3
                                error('Dữ liệu không đủ để tính trung tâm')
                            end
                            if idx(1) > 1 && idx(1) < length(xi)
                                % Sử dụng 3 điểm để tính trung tâm
                                x_minus = xi(idx(1) - 1);
                                x_plus = xi(idx(1) + 1);
                                y_minus = yi(idx(1) - 1);
                                y_plus = yi(idx(1) + 1);
                                result = (y_plus - y_minus)/(x_plus - x_minus);
                            else
                                error('Không thể tính đạo hàm trung tâm với dữ liệu giới hạn!');
                            end
                    end
                elseif ~isempty(func) && ~isnan(h)
                    % Trường hợp nhập hàm số và bước h
                    try
                        f = str2func(['@(x)',func]);
                    catch
                        error('Hàm số không hợp lệ!');
                    end
                    if h<=0
                        error('Bước h phải lớn hơn 0!');
                    end
                 % Tính đạo hàm theo phương pháp
                 switch selectedButton.Text
                     case 'Xấp xỉ tiến'
                         if strcmp(selectedButton2, 'O(h)')
                             result = (f(xdaoham+h) - f(xdaoham))/h;
                         else
                             result = (-f(xdaoham+2*h) + 4*f(xdaoham+h) - 3*f(xdaoham))/(2*h);
                         end
                     case 'Lùi'
                         if strcmp(selectedButton2, 'O(h)')
                            result = (f(xdaoham) - f(xdaoham - h)) / h; % Sai số bậc 1
                        else
                            result = (3*f(xdaoham) - 4*f(xdaoham - h) + f(xdaoham - 2*h)) / (2*h); % Sai số bậc 2
                         end    
                    case 'Trung tâm'
                        if strcmp(selectedButton2, 'O(h)')
                            result = (f(xdaoham + h) - f(xdaoham - h)) / (2*h); % Sai số bậc 1
                        else
                            result = (-f(xdaoham + 2*h) + 8*f(xdaoham + h) - 8*f(xdaoham - h) + f(xdaoham - 2*h)) / (12*h); % Sai số bậc 2
                        end 
                 end
                else
                    error('Hãy nhập đầy đủ dữ liệu x, y hoặc hàm số, bước h và giá trị cần tính đạo hàm!');
                end
                app.KqdaohamEditField.Value = result;
            catch ME
                % Hiển thị thông báo lỗi nếu có
                uialert(app.UIFigure, 'Có lỗi xảy ra! Vui lòng kiểm tra lại dữ liệu.', 'Lỗi');
                disp(ME.message); % Ghi lỗi vào Command Window để dễ dàng kiểm tra
            end
        end

        % Selection changed function: ChontichphanButtonGroup
        function ChontichphanButtonGroupSelectionChanged(app, event)
            selectedButton = app.ChontichphanButtonGroup.SelectedObject;
            try
                xi = str2num(app.NhapxEditField4.Value);
                yi = str2num(app.NhapyEditField4.Value);
                func = app.Nhapham2EditField.Value;
                can = str2num(app.NhapcanEditField.Value);
                if ~isempty(can) && length(can) == 2
                    a = can(1);
                    b = can(2);
                else
                    a = NaN; 
                    b = NaN;
                end
                N = app.NhapNEditField.Value;
                I = 0;
                if ~isempty(xi) && ~isempty(yi)
                    if length(xi) ~= length(yi)
                        error('Dữ liệu x và y phải có cùng số lượng phần tử!');
                    end
                    n = length(xi) - 1;
                    h = xi(2) - xi(1);
                    if any(abs(diff(xi) - h) > 1e-6)
                        error('Các điểm x phải cách đều nhau!');
                    end
                    switch selectedButton.Text
                        case 'Hình thang'
                            I = h*(2*sum(yi) - (yi(1) + yi(end)))/2;
                        case 'Simpson 1/3'
                            if mod(n, 2) ~= 0
                                error('Số đoạn chia phải là số chẵn!')
                            end
                            I = (h/3)*(yi(1) + yi(end) + 4*sum(yi(2:2:end-1)) + 2*sum(yi(3:2:end-2)));
                        case 'Simpson 3/8'
                            if mod(n,3) ~= 0
                                error('Số đoạn chia phải là bội số của 3!');
                            end
                            I = (3/8)*h*(yi(1) + yi(end) + 3*sum(yi(2:3:end)) + 3*sum(yi(3:3:end-1)) + 2*sum(yi(4:3:end-3)));
                    end
                elseif ~isempty(func) && ~isnan(a) && ~isnan(b)
                    try
                        fx = str2func(['@(x)', func]);                                            
                        if N<=0
                            error('N phải lớn hơn 0!');
                        end
                        switch selectedButton.Text
                            case 'Hình thang'
                                I = app.tichphanhinhthang(fx,a,b,N);
                            case 'Simpson 1/3'
                                I = app.tichphanSimpson1_3(fx,a,b,N);
                            case 'Simpson 3/8'
                                I = app.tichphanSimpson3_8(fx,a,b,N);
                        end
                    catch
                        error('Hàm số không hợp lệ! Bạn cần ghi hàm dưới dạng mảng.');
                    end
                else
                    error('Hãy nhập đầy đủ dữ liệu x, y hoặc hàm số, cận và N!');
                end
                app.KqtichphanEditField.Value = I;
            catch ME
                % Hiển thị thông báo lỗi nếu có
                uialert(app.UIFigure, 'Có lỗi xảy ra! Vui lòng kiểm tra lại dữ liệu.', 'Lỗi');
                disp(ME.message); % Ghi lỗi vào Command Window để dễ dàng kiểm tra
            end
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 778 480];
            app.UIFigure.Name = 'MATLAB App';

            % Create TabGroup
            app.TabGroup = uitabgroup(app.UIFigure);
            app.TabGroup.Position = [1 1 779 480];

            % Create NghiemTab
            app.NghiemTab = uitab(app.TabGroup);
            app.NghiemTab.Title = 'Nghiệm';

            % Create UIAxes
            app.UIAxes = uiaxes(app.NghiemTab);
            title(app.UIAxes, 'Đồ thị hàm số')
            xlabel(app.UIAxes, 'x')
            ylabel(app.UIAxes, 'y')
            zlabel(app.UIAxes, 'Z')
            app.UIAxes.XGrid = 'on';
            app.UIAxes.YGrid = 'on';
            app.UIAxes.Position = [324 172 454 236];

            % Create NhpphngtrnhEditFieldLabel
            app.NhpphngtrnhEditFieldLabel = uilabel(app.NghiemTab);
            app.NhpphngtrnhEditFieldLabel.HorizontalAlignment = 'right';
            app.NhpphngtrnhEditFieldLabel.Position = [69 388 107 22];
            app.NhpphngtrnhEditFieldLabel.Text = 'Nhập phương trình';

            % Create PhuongtrinhEditField
            app.PhuongtrinhEditField = uieditfield(app.NghiemTab, 'text');
            app.PhuongtrinhEditField.Position = [191 388 100 22];

            % Create NhpkhongphnlynghimEditFieldLabel
            app.NhpkhongphnlynghimEditFieldLabel = uilabel(app.NghiemTab);
            app.NhpkhongphnlynghimEditFieldLabel.HorizontalAlignment = 'right';
            app.NhpkhongphnlynghimEditFieldLabel.Position = [15 351 161 22];
            app.NhpkhongphnlynghimEditFieldLabel.Text = 'Nhập khoảng phân ly nghiệm';

            % Create PhanlyEditField
            app.PhanlyEditField = uieditfield(app.NghiemTab, 'text');
            app.PhanlyEditField.Position = [191 351 100 22];

            % Create PhuongphapButtonGroup
            app.PhuongphapButtonGroup = uibuttongroup(app.NghiemTab);
            app.PhuongphapButtonGroup.SelectionChangedFcn = createCallbackFcn(app, @PhuongphapButtonGroupSelectionChanged, true);
            app.PhuongphapButtonGroup.Title = 'Chọn phương pháp tìm';
            app.PhuongphapButtonGroup.Position = [110 143 140 153];

            % Create ChiadoiButton
            app.ChiadoiButton = uitogglebutton(app.PhuongphapButtonGroup);
            app.ChiadoiButton.Text = 'Chia đôi';
            app.ChiadoiButton.Position = [11 99 119 23];
            app.ChiadoiButton.Value = true;

            % Create LapButton
            app.LapButton = uitogglebutton(app.PhuongphapButtonGroup);
            app.LapButton.Text = 'Lặp';
            app.LapButton.Position = [11 70 119 23];

            % Create NewtonButton
            app.NewtonButton = uitogglebutton(app.PhuongphapButtonGroup);
            app.NewtonButton.Text = 'Newton (tiếp tuyến)';
            app.NewtonButton.Position = [11 42 119 23];

            % Create DaycungButton
            app.DaycungButton = uitogglebutton(app.PhuongphapButtonGroup);
            app.DaycungButton.Text = 'Dây cung';
            app.DaycungButton.Position = [11 10 119 23];

            % Create KtqunghimEditFieldLabel
            app.KtqunghimEditFieldLabel = uilabel(app.NghiemTab);
            app.KtqunghimEditFieldLabel.HorizontalAlignment = 'right';
            app.KtqunghimEditFieldLabel.Position = [86 95 89 22];
            app.KtqunghimEditFieldLabel.Text = 'Kết quả nghiệm';

            % Create KetquaEditField
            app.KetquaEditField = uieditfield(app.NghiemTab, 'text');
            app.KetquaEditField.Editable = 'off';
            app.KetquaEditField.Position = [190 95 100 22];

            % Create SlnlpEditFieldLabel
            app.SlnlpEditFieldLabel = uilabel(app.NghiemTab);
            app.SlnlpEditFieldLabel.HorizontalAlignment = 'right';
            app.SlnlpEditFieldLabel.Position = [418 95 58 22];
            app.SlnlpEditFieldLabel.Text = 'Số lần lặp';

            % Create SolanEditField
            app.SolanEditField = uieditfield(app.NghiemTab, 'text');
            app.SolanEditField.Editable = 'off';
            app.SolanEditField.Position = [491 95 100 22];

            % Create SaisLabel
            app.SaisLabel = uilabel(app.NghiemTab);
            app.SaisLabel.HorizontalAlignment = 'right';
            app.SaisLabel.Position = [137 315 38 22];
            app.SaisLabel.Text = 'Sai số';

            % Create SaisoEditField
            app.SaisoEditField = uieditfield(app.NghiemTab, 'text');
            app.SaisoEditField.Position = [190 315 100 22];

            % Create NoisuyTab
            app.NoisuyTab = uitab(app.TabGroup);
            app.NoisuyTab.Title = 'Nội suy';

            % Create UIAxes2
            app.UIAxes2 = uiaxes(app.NoisuyTab);
            title(app.UIAxes2, 'Đồ thị hàm số nội suy và dữ liệu thực')
            xlabel(app.UIAxes2, 'X')
            ylabel(app.UIAxes2, 'Y')
            zlabel(app.UIAxes2, 'Z')
            app.UIAxes2.XGrid = 'on';
            app.UIAxes2.YGrid = 'on';
            app.UIAxes2.Position = [298 176 465 232];

            % Create NhpdliuxEditFieldLabel
            app.NhpdliuxEditFieldLabel = uilabel(app.NoisuyTab);
            app.NhpdliuxEditFieldLabel.HorizontalAlignment = 'right';
            app.NhpdliuxEditFieldLabel.Position = [64 388 83 22];
            app.NhpdliuxEditFieldLabel.Text = 'Nhập dữ liệu x';

            % Create NhapxEditField
            app.NhapxEditField = uieditfield(app.NoisuyTab, 'text');
            app.NhapxEditField.Position = [162 388 100 22];

            % Create NhpdliuyEditFieldLabel
            app.NhpdliuyEditFieldLabel = uilabel(app.NoisuyTab);
            app.NhpdliuyEditFieldLabel.HorizontalAlignment = 'right';
            app.NhpdliuyEditFieldLabel.Position = [64 353 83 22];
            app.NhpdliuyEditFieldLabel.Text = 'Nhập dữ liệu y';

            % Create NhapyEditField
            app.NhapyEditField = uieditfield(app.NoisuyTab, 'text');
            app.NhapyEditField.Position = [162 353 100 22];

            % Create NoisuyButtonGroup
            app.NoisuyButtonGroup = uibuttongroup(app.NoisuyTab);
            app.NoisuyButtonGroup.SelectionChangedFcn = createCallbackFcn(app, @NoisuyButtonGroupSelectionChanged, true);
            app.NoisuyButtonGroup.Title = 'Chọn phương pháp nội suy';
            app.NoisuyButtonGroup.Position = [86 235 176 93];

            % Create NewtonButton2
            app.NewtonButton2 = uitogglebutton(app.NoisuyButtonGroup);
            app.NewtonButton2.Text = 'Newton';
            app.NewtonButton2.Position = [38 7 100 23];
            app.NewtonButton2.Value = true;

            % Create LagrangeButton
            app.LagrangeButton = uitogglebutton(app.NoisuyButtonGroup);
            app.LagrangeButton.Text = 'Lagrange';
            app.LagrangeButton.Position = [38 43 100 23];

            % Create NhpgitrcnnisuyEditFieldLabel
            app.NhpgitrcnnisuyEditFieldLabel = uilabel(app.NoisuyTab);
            app.NhpgitrcnnisuyEditFieldLabel.HorizontalAlignment = 'right';
            app.NhpgitrcnnisuyEditFieldLabel.Position = [19 95 130 22];
            app.NhpgitrcnnisuyEditFieldLabel.Text = 'Nhập giá trị cần nội suy';

            % Create NhapnoisuyEditField
            app.NhapnoisuyEditField = uieditfield(app.NoisuyTab, 'numeric');
            app.NhapnoisuyEditField.Position = [164 95 100 22];

            % Create KtqunisuyEditFieldLabel
            app.KtqunisuyEditFieldLabel = uilabel(app.NoisuyTab);
            app.KtqunisuyEditFieldLabel.HorizontalAlignment = 'right';
            app.KtqunisuyEditFieldLabel.Position = [339 95 88 22];
            app.KtqunisuyEditFieldLabel.Text = 'Kết quả nội suy';

            % Create KqnoisuyEditField
            app.KqnoisuyEditField = uieditfield(app.NoisuyTab, 'numeric');
            app.KqnoisuyEditField.Editable = 'off';
            app.KqnoisuyEditField.Position = [442 95 100 22];

            % Create KtquathcnisuyEditFieldLabel
            app.KtquathcnisuyEditFieldLabel = uilabel(app.NoisuyTab);
            app.KtquathcnisuyEditFieldLabel.HorizontalAlignment = 'right';
            app.KtquathcnisuyEditFieldLabel.Position = [16 156 132 22];
            app.KtquathcnisuyEditFieldLabel.Text = 'Kết quả đa thức nội suy';

            % Create KqdathucEditField
            app.KqdathucEditField = uieditfield(app.NoisuyTab, 'text');
            app.KqdathucEditField.Editable = 'off';
            app.KqdathucEditField.Position = [163 143 313 35];

            % Create HoiquyTab
            app.HoiquyTab = uitab(app.TabGroup);
            app.HoiquyTab.Title = 'Hồi quy';

            % Create UIAxes3
            app.UIAxes3 = uiaxes(app.HoiquyTab);
            title(app.UIAxes3, 'Đồ thị hàm hồi quy và dữ liệu thực')
            xlabel(app.UIAxes3, 'x')
            ylabel(app.UIAxes3, 'y')
            zlabel(app.UIAxes3, 'Z')
            app.UIAxes3.XGrid = 'on';
            app.UIAxes3.YGrid = 'on';
            app.UIAxes3.Position = [364 176 400 260];

            % Create NhpdliuxEditFieldLabel_2
            app.NhpdliuxEditFieldLabel_2 = uilabel(app.HoiquyTab);
            app.NhpdliuxEditFieldLabel_2.HorizontalAlignment = 'right';
            app.NhpdliuxEditFieldLabel_2.Position = [93 395 83 22];
            app.NhpdliuxEditFieldLabel_2.Text = 'Nhập dữ liệu x';

            % Create NhapxEditField2
            app.NhapxEditField2 = uieditfield(app.HoiquyTab, 'text');
            app.NhapxEditField2.Position = [191 395 100 22];

            % Create NhpdliuyEditFieldLabel_2
            app.NhpdliuyEditFieldLabel_2 = uilabel(app.HoiquyTab);
            app.NhpdliuyEditFieldLabel_2.HorizontalAlignment = 'right';
            app.NhpdliuyEditFieldLabel_2.Position = [93 360 83 22];
            app.NhpdliuyEditFieldLabel_2.Text = 'Nhập dữ liệu y';

            % Create NhapyEditField2
            app.NhapyEditField2 = uieditfield(app.HoiquyTab, 'text');
            app.NhapyEditField2.Position = [191 360 100 22];

            % Create HoiquyButtonGroup
            app.HoiquyButtonGroup = uibuttongroup(app.HoiquyTab);
            app.HoiquyButtonGroup.SelectionChangedFcn = createCallbackFcn(app, @HoiquyButtonGroupSelectionChanged, true);
            app.HoiquyButtonGroup.Title = 'Chọn phương pháp hồi quy';
            app.HoiquyButtonGroup.Position = [175 208 164 129];

            % Create TuyentinhButton
            app.TuyentinhButton = uitogglebutton(app.HoiquyButtonGroup);
            app.TuyentinhButton.Text = 'Tuyến tính';
            app.TuyentinhButton.Position = [32 74 100 23];
            app.TuyentinhButton.Value = true;

            % Create HammuButton
            app.HammuButton = uitogglebutton(app.HoiquyButtonGroup);
            app.HammuButton.Text = 'Hàm mũ';
            app.HammuButton.Position = [32 43 100 23];

            % Create LogaritButton
            app.LogaritButton = uitogglebutton(app.HoiquyButtonGroup);
            app.LogaritButton.Text = 'Logarit';
            app.LogaritButton.Position = [32 10 100 23];

            % Create KtquphngtrnhhiquyEditFieldLabel
            app.KtquphngtrnhhiquyEditFieldLabel = uilabel(app.HoiquyTab);
            app.KtquphngtrnhhiquyEditFieldLabel.HorizontalAlignment = 'right';
            app.KtquphngtrnhhiquyEditFieldLabel.Position = [19 127 162 22];
            app.KtquphngtrnhhiquyEditFieldLabel.Text = 'Kết quả phương trình hồi quy';

            % Create KqpthoiquyEditField
            app.KqpthoiquyEditField = uieditfield(app.HoiquyTab, 'text');
            app.KqpthoiquyEditField.Editable = 'off';
            app.KqpthoiquyEditField.Position = [191 122 175 32];

            % Create NhpgitrcndonEditFieldLabel
            app.NhpgitrcndonEditFieldLabel = uilabel(app.HoiquyTab);
            app.NhpgitrcndonEditFieldLabel.HorizontalAlignment = 'right';
            app.NhpgitrcndonEditFieldLabel.Position = [38 83 137 22];
            app.NhpgitrcndonEditFieldLabel.Text = 'Nhập giá trị cần dự đoán';

            % Create NhapdudoanEditField
            app.NhapdudoanEditField = uieditfield(app.HoiquyTab, 'numeric');
            app.NhapdudoanEditField.Position = [190 83 100 22];

            % Create KtqudonEditFieldLabel
            app.KtqudonEditFieldLabel = uilabel(app.HoiquyTab);
            app.KtqudonEditFieldLabel.HorizontalAlignment = 'right';
            app.KtqudonEditFieldLabel.Position = [440 83 94 22];
            app.KtqudonEditFieldLabel.Text = 'Kết quả dự đoán';

            % Create KqdudoanEditField
            app.KqdudoanEditField = uieditfield(app.HoiquyTab, 'numeric');
            app.KqdudoanEditField.Editable = 'off';
            app.KqdudoanEditField.Position = [549 83 100 22];

            % Create DaohamTab
            app.DaohamTab = uitab(app.TabGroup);
            app.DaohamTab.Title = 'Đạo hàm';

            % Create NhpdliuxEditFieldLabel_3
            app.NhpdliuxEditFieldLabel_3 = uilabel(app.DaohamTab);
            app.NhpdliuxEditFieldLabel_3.HorizontalAlignment = 'right';
            app.NhpdliuxEditFieldLabel_3.Position = [79 375 83 22];
            app.NhpdliuxEditFieldLabel_3.Text = 'Nhập dữ liệu x';

            % Create NhapxEditField3
            app.NhapxEditField3 = uieditfield(app.DaohamTab, 'text');
            app.NhapxEditField3.Position = [177 375 100 22];

            % Create NhpdliuyEditFieldLabel_3
            app.NhpdliuyEditFieldLabel_3 = uilabel(app.DaohamTab);
            app.NhpdliuyEditFieldLabel_3.HorizontalAlignment = 'right';
            app.NhpdliuyEditFieldLabel_3.Position = [79 340 83 22];
            app.NhpdliuyEditFieldLabel_3.Text = 'Nhập dữ liệu y';

            % Create NhapyEditField3
            app.NhapyEditField3 = uieditfield(app.DaohamTab, 'text');
            app.NhapyEditField3.Position = [177 340 100 22];

            % Create HocLabel
            app.HocLabel = uilabel(app.DaohamTab);
            app.HocLabel.FontSize = 18;
            app.HocLabel.Position = [150 290 47 32];
            app.HocLabel.Text = 'Hoặc';

            % Create NhphmsEditFieldLabel
            app.NhphmsEditFieldLabel = uilabel(app.DaohamTab);
            app.NhphmsEditFieldLabel.HorizontalAlignment = 'right';
            app.NhphmsEditFieldLabel.Position = [87 246 76 22];
            app.NhphmsEditFieldLabel.Text = 'Nhập hàm số';

            % Create NhaphamEditField
            app.NhaphamEditField = uieditfield(app.DaohamTab, 'text');
            app.NhaphamEditField.Position = [177 242 188 30];

            % Create ChonsaisoButtonGroup
            app.ChonsaisoButtonGroup = uibuttongroup(app.DaohamTab);
            app.ChonsaisoButtonGroup.SelectionChangedFcn = createCallbackFcn(app, @ChondaohamButtonGroupSelectionChanged, true);
            app.ChonsaisoButtonGroup.Title = 'Chọn giá trị sai sô';
            app.ChonsaisoButtonGroup.Position = [117 80 123 105];

            % Create OhButton
            app.OhButton = uitogglebutton(app.ChonsaisoButtonGroup);
            app.OhButton.Text = 'O(h)';
            app.OhButton.Position = [11 51 100 23];
            app.OhButton.Value = true;

            % Create Oh2Button
            app.Oh2Button = uitogglebutton(app.ChonsaisoButtonGroup);
            app.Oh2Button.Text = 'O(h^2)';
            app.Oh2Button.Position = [11 16 100 23];

            % Create ChondaohamButtonGroup
            app.ChondaohamButtonGroup = uibuttongroup(app.DaohamTab);
            app.ChondaohamButtonGroup.SelectionChangedFcn = createCallbackFcn(app, @ChondaohamButtonGroupSelectionChanged, true);
            app.ChondaohamButtonGroup.Title = 'Chọn phương pháp đạo hàm';
            app.ChondaohamButtonGroup.Position = [485 258 185 128];

            % Create XpxtinButton
            app.XpxtinButton = uitogglebutton(app.ChondaohamButtonGroup);
            app.XpxtinButton.Text = 'Xấp xỉ tiến';
            app.XpxtinButton.Position = [43 77 100 23];
            app.XpxtinButton.Value = true;

            % Create LiButton
            app.LiButton = uitogglebutton(app.ChondaohamButtonGroup);
            app.LiButton.Text = 'Lùi';
            app.LiButton.Position = [43 43 100 23];

            % Create TrungtmButton
            app.TrungtmButton = uitogglebutton(app.ChondaohamButtonGroup);
            app.TrungtmButton.Text = 'Trung tâm';
            app.TrungtmButton.Position = [43 10 100 23];

            % Create KtquEditFieldLabel
            app.KtquEditFieldLabel = uilabel(app.DaohamTab);
            app.KtquEditFieldLabel.HorizontalAlignment = 'right';
            app.KtquEditFieldLabel.Position = [528 150 46 22];
            app.KtquEditFieldLabel.Text = 'Kết quả';

            % Create KqdaohamEditField
            app.KqdaohamEditField = uieditfield(app.DaohamTab, 'numeric');
            app.KqdaohamEditField.Editable = 'off';
            app.KqdaohamEditField.Position = [589 150 100 22];

            % Create NhpgitrcntnhohmEditField_2Label
            app.NhpgitrcntnhohmEditField_2Label = uilabel(app.DaohamTab);
            app.NhpgitrcntnhohmEditField_2Label.HorizontalAlignment = 'right';
            app.NhpgitrcntnhohmEditField_2Label.Position = [412 205 162 22];
            app.NhpgitrcntnhohmEditField_2Label.Text = 'Nhập giá trị cần tính đạo hàm';

            % Create NhapdaohamEditField
            app.NhapdaohamEditField = uieditfield(app.DaohamTab, 'numeric');
            app.NhapdaohamEditField.Position = [589 205 100 22];

            % Create NhpbchEditFieldLabel
            app.NhpbchEditFieldLabel = uilabel(app.DaohamTab);
            app.NhpbchEditFieldLabel.HorizontalAlignment = 'right';
            app.NhpbchEditFieldLabel.Position = [86 208 76 22];
            app.NhpbchEditFieldLabel.Text = 'Nhập bước h';

            % Create NhaphEditField
            app.NhaphEditField = uieditfield(app.DaohamTab, 'numeric');
            app.NhaphEditField.Position = [177 208 100 22];

            % Create TichphanTab
            app.TichphanTab = uitab(app.TabGroup);
            app.TichphanTab.Title = 'Tích phân';

            % Create NhpdliuxEditFieldLabel_4
            app.NhpdliuxEditFieldLabel_4 = uilabel(app.TichphanTab);
            app.NhpdliuxEditFieldLabel_4.HorizontalAlignment = 'right';
            app.NhpdliuxEditFieldLabel_4.Position = [73 376 83 22];
            app.NhpdliuxEditFieldLabel_4.Text = 'Nhập dữ liệu x';

            % Create NhapxEditField4
            app.NhapxEditField4 = uieditfield(app.TichphanTab, 'text');
            app.NhapxEditField4.Position = [171 376 100 22];

            % Create NhpdliuyEditFieldLabel_4
            app.NhpdliuyEditFieldLabel_4 = uilabel(app.TichphanTab);
            app.NhpdliuyEditFieldLabel_4.HorizontalAlignment = 'right';
            app.NhpdliuyEditFieldLabel_4.Position = [73 341 83 22];
            app.NhpdliuyEditFieldLabel_4.Text = 'Nhập dữ liệu y';

            % Create NhapyEditField4
            app.NhapyEditField4 = uieditfield(app.TichphanTab, 'text');
            app.NhapyEditField4.Position = [171 341 100 22];

            % Create HocLabel_2
            app.HocLabel_2 = uilabel(app.TichphanTab);
            app.HocLabel_2.FontSize = 18;
            app.HocLabel_2.Position = [148 284 47 32];
            app.HocLabel_2.Text = 'Hoặc';

            % Create NhphmsEditFieldLabel_2
            app.NhphmsEditFieldLabel_2 = uilabel(app.TichphanTab);
            app.NhphmsEditFieldLabel_2.HorizontalAlignment = 'right';
            app.NhphmsEditFieldLabel_2.Position = [47 237 76 22];
            app.NhphmsEditFieldLabel_2.Text = 'Nhập hàm số';

            % Create Nhapham2EditField
            app.Nhapham2EditField = uieditfield(app.TichphanTab, 'text');
            app.Nhapham2EditField.Position = [137 233 242 30];

            % Create ChontichphanButtonGroup
            app.ChontichphanButtonGroup = uibuttongroup(app.TichphanTab);
            app.ChontichphanButtonGroup.SelectionChangedFcn = createCallbackFcn(app, @ChontichphanButtonGroupSelectionChanged, true);
            app.ChontichphanButtonGroup.Title = 'Chọn phương pháp Tích phân';
            app.ChontichphanButtonGroup.Position = [514 289 171 136];

            % Create HnhthangButton
            app.HnhthangButton = uitogglebutton(app.ChontichphanButtonGroup);
            app.HnhthangButton.Text = 'Hình thang';
            app.HnhthangButton.Position = [36 79 100 23];
            app.HnhthangButton.Value = true;

            % Create Simpson13Button
            app.Simpson13Button = uitogglebutton(app.ChontichphanButtonGroup);
            app.Simpson13Button.Text = 'Simpson 1/3';
            app.Simpson13Button.Position = [36 46 100 23];

            % Create Simpson38Button
            app.Simpson38Button = uitogglebutton(app.ChontichphanButtonGroup);
            app.Simpson38Button.Text = 'Simpson 3/8';
            app.Simpson38Button.Position = [36 9 100 23];

            % Create NhpNEditFieldLabel
            app.NhpNEditFieldLabel = uilabel(app.TichphanTab);
            app.NhpNEditFieldLabel.HorizontalAlignment = 'right';
            app.NhpNEditFieldLabel.Position = [543 205 46 22];
            app.NhpNEditFieldLabel.Text = 'Nhập N';

            % Create NhapNEditField
            app.NhapNEditField = uieditfield(app.TichphanTab, 'numeric');
            app.NhapNEditField.Position = [604 205 100 22];

            % Create KtquEditFieldLabel_2
            app.KtquEditFieldLabel_2 = uilabel(app.TichphanTab);
            app.KtquEditFieldLabel_2.HorizontalAlignment = 'right';
            app.KtquEditFieldLabel_2.Position = [284 132 46 22];
            app.KtquEditFieldLabel_2.Text = 'Kết quả';

            % Create KqtichphanEditField
            app.KqtichphanEditField = uieditfield(app.TichphanTab, 'numeric');
            app.KqtichphanEditField.Editable = 'off';
            app.KqtichphanEditField.Position = [345 132 179 22];

            % Create NhpcntnhtchphnEditFieldLabel
            app.NhpcntnhtchphnEditFieldLabel = uilabel(app.TichphanTab);
            app.NhpcntnhtchphnEditFieldLabel.HorizontalAlignment = 'right';
            app.NhpcntnhtchphnEditFieldLabel.Position = [440 243 149 22];
            app.NhpcntnhtchphnEditFieldLabel.Text = 'Nhập cận để tính tích phân';

            % Create NhapcanEditField
            app.NhapcanEditField = uieditfield(app.TichphanTab, 'text');
            app.NhapcanEditField.Position = [604 243 100 22];

            % Create GioithieunhomTab
            app.GioithieunhomTab = uitab(app.TabGroup);
            app.GioithieunhomTab.Title = 'Giới thiệu Nhóm';

            % Create Nhm7Label
            app.Nhm7Label = uilabel(app.GioithieunhomTab);
            app.Nhm7Label.HorizontalAlignment = 'center';
            app.Nhm7Label.FontSize = 18;
            app.Nhm7Label.FontWeight = 'bold';
            app.Nhm7Label.Position = [337 402 71 23];
            app.Nhm7Label.Text = 'Nhóm 7';

            % Create ThnhvinLabel
            app.ThnhvinLabel = uilabel(app.GioithieunhomTab);
            app.ThnhvinLabel.FontWeight = 'bold';
            app.ThnhvinLabel.Position = [64 351 76 22];
            app.ThnhvinLabel.Text = 'Thành viên: ';

            % Create NguynQucTrng22200172Label
            app.NguynQucTrng22200172Label = uilabel(app.GioithieunhomTab);
            app.NguynQucTrng22200172Label.Position = [146 351 210 22];
            app.NguynQucTrng22200172Label.Text = 'Nguyễn Quốc Trường - 22200172';

            % Create NhimvLabel
            app.NhimvLabel = uilabel(app.GioithieunhomTab);
            app.NhimvLabel.HorizontalAlignment = 'center';
            app.NhimvLabel.FontWeight = 'bold';
            app.NhimvLabel.Position = [68 321 63 22];
            app.NhimvLabel.Text = 'Nhiệm vụ:';

            % Create ThitkvcodetonbnthchnhphngphptnhLabel
            app.ThitkvcodetonbnthchnhphngphptnhLabel = uilabel(app.GioithieunhomTab);
            app.ThitkvcodetonbnthchnhphngphptnhLabel.Position = [146 321 332 22];
            app.ThitkvcodetonbnthchnhphngphptnhLabel.Text = 'Thiết kế và code toàn bộ đồ án thực hành phương pháp tính';

            % Create GiithiuLabel
            app.GiithiuLabel = uilabel(app.GioithieunhomTab);
            app.GiithiuLabel.HorizontalAlignment = 'center';
            app.GiithiuLabel.FontWeight = 'bold';
            app.GiithiuLabel.Position = [67 295 66 22];
            app.GiithiuLabel.Text = 'Giới thiệu:';

            % Create nthitkapptrnMatlabvccphntnhtonhcnhNghimNisuyHiquyohmTchphnLabel
            app.nthitkapptrnMatlabvccphntnhtonhcnhNghimNisuyHiquyohmTchphnLabel = uilabel(app.GioithieunhomTab);
            app.nthitkapptrnMatlabvccphntnhtonhcnhNghimNisuyHiquyohmTchphnLabel.Position = [146 295 620 22];
            app.nthitkapptrnMatlabvccphntnhtonhcnhNghimNisuyHiquyohmTchphnLabel.Text = 'Đồ án thiết kế app trên Matlab về các phần tính toán đã học như: Nghiệm, Nội suy, Hồi quy, Đạo hàm, Tích phân.';

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = project_exported

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end