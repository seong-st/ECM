classdef CellEcm
    properties
        current
        time
        voltage
        id1
        id2
        id3
        id4
        eta_chg
        eta_dis
        q_cell
    end
    
    methods
        function obj = CellEcm(data, params)
            % Initialize with HPPC data and model parameters
            obj.current = data.current;
            obj.time = data.time;
            obj.voltage = data.voltage;
            [obj.id1, obj.id2, obj.id3, obj.id4] = data.get_indices();

            obj.eta_chg = params.eta_chg;
            obj.eta_dis = params.eta_dis;
            obj.q_cell = params.q_cell;

%             obj.time = obj.time(1:10:end);
%             obj.voltage = obj.voltage(1:10:end);
%             obj.current = obj.current(1:10:end);
        end
        
        function soc = soc(obj)
            current = obj.current;
            time = obj.time;
            q = obj.q_cell * 3600;
            dt = diff(time);

            nc = length(current);
            soc = ones(nc, 1);

            for k = 2:nc
                i = current(k);
                if i < 0
                    eta = obj.eta_chg;
                else
                    eta = obj.eta_dis;
                end
                soc(k) = soc(k - 1) - ((eta * i * dt(k - 1)) / q);
            end
%             disp('Calculated SOC:');
%             disp(soc); %
        end
        
        function soc = soc2(obj, data)
            current = data.current;
            time = data.time;
            q = obj.q_cell * 3600;
            dt = diff(time);

            nc = length(current);
            soc = ones(nc, 1);

            for k = 2:nc
                i = current(k);
                if i < 0
                    eta = obj.eta_chg;
                else
                    eta = obj.eta_dis;
                end
                soc(k) = soc(k - 1) - ((eta * i * dt(k - 1)) / q);
            end
        end
        
        function [ocv, i_pts, t_pts, v_pts, z_pts] = ocv(obj, soc, pts, vz_pts)
            if nargin < 3
                pts = false;
            end
            if nargin < 4
                vz_pts = [];
            end
            
            nc = length(obj.id1);
            v_pts = obj.voltage(1);
            z_pts = soc(1);
            i_pts = obj.current(1);
            t_pts = obj.time(1);

            v_array = obj.voltage;
            current_array = obj.current;
            time_array = obj.time;
             
%             disp('Initial values:');
%             disp(['v_pts: ', num2str(v_pts)]);
%             disp(['z_pts: ', num2str(z_pts)]);

            if pts 
%                 v_pts = [];
%                 z_pts = [];
%                 i_pts = [];
%                 t_pts = [];
                v_pts = [v_array(1)];
                z_pts = [soc(1)];
                i_pts = [current_array(1)];
                t_pts = [time_array(1)];

                for k = 2:nc
                    aa = obj.id1(k);
                    v_pts = [v_pts; v_array(aa)];
                    z_pts = [z_pts; soc(aa)];
                    i_pts = [i_pts; current_array(aa)];
                    t_pts = [t_pts; time_array(aa)];                
                    
                    disp([num2str(v_pts)]);
                end
                

                v_pts = [v_pts; v_array(end)];
                z_pts = [z_pts; soc(end)];
                i_pts = [i_pts; current_array(end)];
                t_pts = [t_pts; time_array(end)];

                [z_pts, unique_idx] = unique(z_pts);
                v_pts = v_pts(unique_idx);

                disp('After loop:');
                disp(['v_pts: ', num2str(v_pts')]);
                disp(['z_pts: ', num2str(z_pts')]);
            

                ocv = interp1(z_pts(end:-1:1), v_pts(end:-1:1), soc);
            elseif ~isempty(vz_pts)
                v_pts = vz_pts{1};
                z_pts = vz_pts{2};

                [z_pts, unique_idx] = unique(z_pts);
                v_pts = v_pts(unique_idx);

                ocv = interp1(z_pts(end:-1:1), v_pts(end:-1:1), soc);
            end

            v_pts = v_pts(:);
            z_pts = z_pts(:);
            disp('Returning values from function:');
            disp(['v_pts: ', num2str(v_pts')]);
            disp(['z_pts: ', num2str(z_pts')]);           

        end
        
        function coeffs = curve_fit_coeff(obj, func, ncoeff)
            rest_start = obj.id4;
            rest_end = obj.id1;
            nrow = length(obj.id4);
            coeffs = zeros(nrow, ncoeff);

            for i = 1:nrow  
                start = rest_start(i);
                end_ = rest_end(i+1);
                t_curve = obj.time(start:end_);
                v_curve = obj.voltage(start:end_);
                
                t_scale = t_curve - t_curve(1);  
                guess = [v_curve(end), 0.002, 0.01, 0.001, 0.01];

                
                t_scale = t_scale(:);
                v_curve = v_curve(:);

                % nlinfit 함수 호출
                popt = nlinfit(t_scale, v_curve, @(b, t) func(t, b(1), b(2), b(3), b(4), b(5)), guess);

               
                if length(popt) == ncoeff
                    coeffs(i, :) = popt;
                else
                    error('nlinfit 결과의 요소 개수가 예상과 일치하지 않습니다.');
                end
            end
        end
        
        function rctau = rctau_ttc(obj, coeffs)
            s1 = obj.id1;  % 휴식 끝
            s2 = obj.id2;  % 펄스 시작
            s3 = obj.id3;  % 펄스 끝
            s4 = obj.id4;  % 휴식 시작
            
            nrow = length(s4);
            num_coeffs = size(coeffs, 1);
            rctau = zeros(num_coeffs, 7);  
                
            disp('coeffs 배열의 크기:');
            disp(size(coeffs));
            disp('coeffs 배열의 내용:');
            disp(coeffs);

            for k = 1:num_coeffs  
                di = abs(obj.current(s1(k)) - obj.current(s2(k)));
                dt = obj.time(s1(k+1)) - obj.time(s4(k));
                dv = abs(obj.voltage(s1(k)) - obj.voltage(s2(k)));
                
                % coeffs의 길이 확인 및 오류 처리
                if size(coeffs, 2) == 5
                    b = coeffs(k, 2);
                    c = coeffs(k, 3);
                    alpha = coeffs(k, 4);
                    beta = coeffs(k, 5);
                else
                    error('coeffs 배열의 요소 개수가 예상과 일치하지 않습니다.');
                end

                tau1 = 1 / alpha;
                tau2 = 1 / beta;
                r0 = dv / di;
                r1 = b / ((1 - exp(-dt / tau1)) * di);
                r2 = c / ((1 - exp(-dt / tau2)) * di);
                c1 = tau1 / r1;
                c2 = tau2 / r2;

                rctau(k, :) = [tau1, tau2, r0, r1, r2, c1, c2];
            end
        end
        function vt = vt(obj, soc, ocv, rctau)
            dt = diff(obj.time);  
            nc = length(obj.current); 
            v0 = zeros(nc, 1);  % v0 배열 초기화
            v1 = zeros(nc, 1);  % v1 배열 초기화
            v2 = zeros(nc, 1);  % v2 배열 초기화

            for k = 2:nc  
                i = obj.current(k);

                % soc에서 파라미터 얻기
                [tau1, tau2, r0, r1, r2] = obj.get_rtau(rctau, soc(k));

                % r0 저항에서의 전압
                v0(k) = r0 * i;

                % c1 커패시터에서의 전압
                tm1 = v1(k - 1) * exp(-dt(k - 1) / tau1);
                tm2 = r1 * (1 - exp(-dt(k - 1) / tau1)) * i;
                v1(k) = tm1 + tm2;

                % c2 커패시터에서의 전압
                tm3 = v2(k - 1) * exp(-dt(k - 1) / tau2);
                tm4 = r2 * (1 - exp(-dt(k - 1) / tau2)) * i;
                v2(k) = tm3 + tm4;
            end

            vt = ocv - v0 - v1 - v2;
        end
        
        function vt2 = vt2(obj, data, soc, ocv, rctau)
            dt = diff(data.time);  
            nc = length(data.current);  
            v0 = zeros(nc, 1);  % v0 배열 초기화
            v1 = zeros(nc, 1);  % v1 배열 초기화
            v2 = zeros(nc, 1);  % v2 배열 초기화
            ocv_new = zeros(nc, 1);
            x = soc;
            y = ocv;
            
            for k = 2:nc  
                i = data.current(k);

                % soc에서 파라미터 얻기
                [tau1, tau2, r0, r1, r2] = obj.get_rtau(rctau, soc(k));

                % r0 저항에서의 전압
                v0(k) = r0 * i;

                % c1 커패시터에서의 전압
                tm1 = v1(k - 1) * exp(-dt(k - 1) / tau1);
                tm2 = r1 * (1 - exp(-dt(k - 1) / tau1)) * i;
                v1(k) = tm1 + tm2;

                % c2 커패시터에서의 전압
                tm3 = v2(k - 1) * exp(-dt(k - 1) / tau2);
                tm4 = r2 * (1 - exp(-dt(k - 1) / tau2)) * i;
                v2(k) = tm3 + tm4;
                
                ocv_new(k) = interp1(x, y, soc(k), 'linear');
                ocv_new(1) = ocv_new(2);
            end

            vt2 = ocv_new - v0 - v1 - v2;
        end

    end
    
    methods (Static)
        function value = func_ttc(t, a, b, c, alpha, beta)
            value = a - b * exp(-alpha * t) - c * exp(-beta * t);
        end
        
        function [tau1, tau2, r0, r1, r2] = get_rtau(rctau, z)
            soc = fliplr(0.0:0.05:1.0);
            [~, idx] = min(abs(soc - z));
            
            idx = min(max(idx, 1), size(rctau, 1));

            tau1 = rctau(idx, 1);
            tau2 = rctau(idx, 2);
            r0 = rctau(idx, 3);
            r1 = rctau(idx, 4);
            r2 = rctau(idx, 5);
        end
    end
end