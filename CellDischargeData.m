classdef CellDischargeData

    % battery discharge data %
    
    properties
        time
        voltage
        current
    end
    
    methods
        function obj = CellDischargeData(path)
            % 생성자: 파일에서 데이터를 읽어와 객체의 속성을 초기화
            data = readtable(path, 'Delimiter', '\t');
            obj.time = data{:, 1};
            obj.voltage = data{:, 2};
            obj.current = data{:, 3};
        end
    end
    
    methods (Static)
        function data = process(path)
            % 원본 방전 데이터의 한 구간을 처리하는 정적 메서드
            data = CellDischargeData(path);
        end
        
        function data = process_discharge_only(path)
            % 방전 데이터의 방전 부분만 처리하는 정적 메서드
            data = CellDischargeData(path);
        end
    end
end

