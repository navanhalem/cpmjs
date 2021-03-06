import resolve from 'rollup-plugin-node-resolve';
import commonjs from 'rollup-plugin-commonjs';
import { eslint } from 'rollup-plugin-eslint';

export default [{
  input: 'app/index.js',
  output: {
    file: 'build/cpm.js',
    format: 'iife',
	name: 'CPM'
  },
  plugins: [
    eslint(),
    resolve(),
    commonjs()
  ]
},{
  input: 'app/index.js',
  output: {
    file: 'build/cpm-cjs.js',
    format: 'cjs',
	name: 'CPM'
  },
  plugins:[
  ]
}];
