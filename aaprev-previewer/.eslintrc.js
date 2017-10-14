module.exports = {
  'env': {
    'browser': true,
    'es6': true
  },
  'extends': [
    'standard',
    'eslint-config-react-app'
  ],
  'parserOptions': {
    'ecmaVersion': 8,
    'ecmaFeatures': {
      'experimentalObjectRestSpread': true,
      'generators': true,
      'jsx': true
    },
    'sourceType': 'module'
  },
  'rules': {
    'indent': [
      'error',
      2,
      {
        'SwitchCase': 1,
        'ignoredNodes': ['JSXAttribute', 'JSXSpreadAttribute']
      }
    ],
    'max-len': [
      'warn',
      100,
      2
    ],
    'react/jsx-indent-props': [2, 1],
    'react/jsx-indent': [2, 2]
  }
}
