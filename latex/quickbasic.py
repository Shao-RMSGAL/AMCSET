from pygments.lexer import RegexLexer, words
from pygments.token import *


class QuickBASICLexer(RegexLexer):
    name = 'QuickBASIC'
    aliases = ['quickbasic']
    filenames = ['*.bas']

    tokens = {
        'root': [
            (r'\b(DEF|DIM|DO|ELSE|END(?:-IF|DO|FUNC|SELECT|SUB|TRY|TYPE|WHILE)?|EXIT(?:-DO|-(?:FOR|FUNCTION|SUB|TRY|WHILE))?|FOR|FUNCTION|IF|LOOP|NEXT|PRINT|RETURN|SELECT|SUB|THEN|TO|UNTIL|WEND|WHILE)\b',Keyword),
            (r'\b(AND|AS|BOX|CASE|CHAIN|CHDIR|CLEAR|CLS|COMMON|CONST|DATA|DECLARE|DEF FN|DECLARE SUB|DEF SEG|DIM SHARED|DRAW|EDIT|ENVIRON|ERASE|FIELD|FILES|FLASH|FONT|GET|INKEY|INPUT|KEY|KILL|LPRINT|LSET|MERGE|MKDIR|NAME|ON ERROR GOTO|ON|OPEN|OPTION|OUT|PAINT|PALETTE|PCOPY|POKE|PRESET|PUT|RANDOMIZE|READ|REDIM|RESET|RESTORE|RESUME|RESUME NEXT|RSET|SAVE|SCREEN|SHELL|SLEEP|SOUND|SPC|SWAP|SYSTEM|TYPE|VIEW|WAIT|WIDTH|WINDOW|WRITE)\b', Keyword.Reserved),
            (r'\b(ARRAY|CALL|CASE|CONST|DO|ELSE|ELSEIF|END|ENDIF|ENDSELECT|ENDSUB|ENDTYPE|ENDWHILE|EXIT|FOR|FUNCTION|IF|LOOP|NEXT|REDO|SELECT|SUB|THEN|TO|UNTIL|WEND|WHILE)\b', Keyword.Pseudo),
            (r'\b[A-Za-z][A-Za-z0-9]*[$%]?\b', Name.Variable),
            (r'\b[A-Za-z_]\w*\b', Name),
            (r'\b\d+\b', Number.Integer),
            (r'\b[\d.]+\s*(?:e[+-]?\d+)?\b', Number.Float),
            (r'"[^"]*"', String),
            (r"'['^\n]*", Comment),
            (r'[()\[\]{}:.,^=]', Punctuation),
            (r'\s+', Text),
            (r'[-+*/!#<>=;]', Operator), #  Include these characters as operators
        ]
    }

# from pygments.lexer import RegexLexer, words
# from pygments.token import *

# class QuickBASICLexer(RegexLexer):
#     name = 'QuickBASIC'
#     aliases = ['quickbasic']
#     filenames = ['*.bas']

#     tokens = {
#         'root': [
#             (r'\b(DEF|DIM|DO|ELSE|END(?:-IF|DO|FUNC|SELECT|SUB|TRY|TYPE|WHILE)?|EXIT(?:-DO|-(?:FOR|FUNCTION|SUB|TRY|WHILE))?|FOR|FUNCTION|IF|LOOP|NEXT|PRINT|RETURN|SELECT|SUB|THEN|TO|UNTIL|WEND|WHILE)\b', Keyword),
#             (r'\b(AND|AS|BOX|CASE|CHAIN|CHDIR|CLEAR|CLS|COMMON|CONST|DATA|DECLARE|DEF FN|DECLARE SUB|DEF SEG|DIM SHARED|DRAW|DRAW|EDIT|ENVIRON|ERASE|FIELD|FILES|FLASH|FONT|GET|INKEY|INPUT|KEY|KILL|LPRINT|LSET|MERGE|MKDIR|NAME|ON ERROR GOTO|ON|OPEN|OPTION|OUT|PAINT|PALETTE|PCOPY|POKE|PRESET|PUT|RANDOMIZE|READ|REDIM|RESET|RESTORE|RESUME|RESUME NEXT|RSET|SAVE|SCREEN|SHELL|SLEEP|SOUND|SPC|SWAP|SYSTEM|TYPE|VIEW|WAIT|WIDTH|WINDOW|WRITE)\b', Keyword.Reserved),
#             (r'\b(ARRAY|CALL|CASE|CONST|DO|ELSE|ELSEIF|END|ENDIF|ENDSELECT|ENDSUB|ENDTYPE|ENDWHILE|EXIT|FOR|FUNCTION|IF|LOOP|NEXT|REDO|SELECT|SUB|THEN|TO|UNTIL|WEND|WHILE)\b', Keyword.Pseudo),
#             (r'\b[A-Za-z][A-Za-z0-9]*[$%]?\b', Name.Variable),
#             (r'\b[A-Za-z_]\w*\b', Name),
#             (r'\b\d+\b', Number.Integer),
#             (r'\b[\d.]+\s*(?:e[+-]?\d+)?\b', Number.Float),
#             (r'"[^"]*"', String),
#             (r"'[^\n]*", Comment),
#             (r'[()\[\]{}:.,^=]', Punctuation),
#             (r'[-+*/!#<>=;]', Operator), #  Include these characters as operators
#             (r'\s+', Text),
#         ]
#     }
