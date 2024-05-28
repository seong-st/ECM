classdef CellHppcData
    properties
        time
        voltage
        current
        flags
    end
    
    methods
        function obj = CellHppcData(path)
            % 생성자: 파일에서 데이터를 읽어와 객체의 속성을 초기화
            data = readtable(path, 'Delimiter', '\t');
            obj.time = data{:, 1};
            obj.voltage = data{:, 2};
            obj.current = data{:, 3};
            obj.flags = data{:, 4};
        end
        
        function [id_end_rest, id_start_pulse, id_end_pulse, id_start_rest] = get_indices(obj)
            % 플래그에 따라 인덱스를 찾는 메서드
            id_end_rest = find(obj.flags == 1);
            id_start_pulse = find(obj.flags == 2);
            id_end_pulse = find(obj.flags == 3);
            id_start_rest = find(obj.flags == 4);
        end
    end
end