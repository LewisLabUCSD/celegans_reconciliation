function num_atoms = create_atom_numbers(metFormula,element,disp_flag)
% num_atoms = create_atom_numbers(metFormula,element,disp_flag)
% gives the number of atoms within a formula
% 
% INPUT:
% metFormula: metabolite formula
% element: which element atom numbers are desired
% disp_flag: display output
% 
% OUTPUT:
% num_atoms: number of atoms
%
% Comment: Accounts for X(Y)n, but nothing more complex..
% Code still needs optimizing

if nargin < 3
    disp_flag = false;
end
all_elements = {'C';'Ca';'Cl';'Co';'Cu';'Fe';'H';'I';'K';'Mg';'Mo';'N';'Na';'O';'P';'R';'S';'Se';'Zn'}';
metFormula_ori = metFormula;
index = regexp(metFormula,element);
index_ = index;
number_of_atoms = 0;
for i=1:length(index_)
    index = index_(i);
    check_elements = setdiff(all_elements,element);
    check_elements = strjoin(check_elements,'|');
    other_index = regexp(metFormula,check_elements);
    index(index==intersect(other_index,index)) = [];
    if ~isempty(index)
%         if i~=1
%             index = index - index_(i-1);
%         end
        if length(index)==1
            if index==length(metFormula) || strcmp(metFormula,element) || ~isempty(regexp(metFormula(index+length(element)),'[A-Z]'))
                num_atoms = num2str(1);
            else
                index2 = min(regexp(metFormula(index+length(element):end),'[A-Z]'));
                if ~isempty(index2)
                    metFormula1 = metFormula(index+length(element):end);
                    num_atoms = metFormula1(1:index2-1);
                else
                    num_atoms = metFormula(index+length(element):end);
                end
            end
        else
            if disp_flag
                fprintf('%s --> Formula not represented correclty.\n',metFormula_ori);
            end
        end
    else
        num_atoms = num2str(0);
    end
    number_of_atoms = number_of_atoms + str2double(num_atoms);
end
num_atoms = num2str(number_of_atoms);